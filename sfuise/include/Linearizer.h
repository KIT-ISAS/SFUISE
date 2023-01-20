#pragma once

#include "utils/common_utils.h"
#include "utils/math_tools.h"
#include "Accumulator.h"
#include "SplineState.h"
#include "Residuals.h"

struct Parameters {

    bool if_opt_g;
    bool if_opt_transform;
    bool if_reject_uwb;
    bool if_reject_uwb_in_optimization;

    double w_uwb;
    double w_acc;
    double w_gyro;
    double w_bias_acc;
    double w_bias_gyro;
    double reject_uwb_thresh;
    double reject_uwb_window_width;

    int control_point_fps;

    Eigen::Vector3d accel_var_inv, gyro_var_inv;
    Eigen::Vector3d bias_accel_var_inv, bias_gyro_var_inv;
    Eigen::Vector3d pos_var_inv;
    Eigen::Vector3d t_nav_uwb_init;
    Eigen::Quaterniond q_nav_uwb_init;

    Eigen::aligned_map<uint16_t, Eigen::Vector3d> anchor_list;
    Eigen::aligned_map<uint16_t, double> toa_offset;

    Parameters() : toa_offset{{(uint16_t)0, 0}} {}

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct Linearizer
{
    static const int POSE_SIZE = 6;
    static const int POS_SIZE = 3;
    static const int POS_OFFSET = 0;
    static const int ROT_SIZE = 3;
    static const int ROT_OFFSET = 3;
    static const int ACCEL_BIAS_SIZE = 3;
    static const int GYRO_BIAS_SIZE = 3;
    static const int BIAS_SIZE = ACCEL_BIAS_SIZE + GYRO_BIAS_SIZE;
    static const int G_SIZE = 2;
    static const int OFFSET_SIZE = 3;
    static const int ROTATION_SIZE = 3;
    static const int BIAS_BLOCK_SIZE = ACCEL_BIAS_SIZE + GYRO_BIAS_SIZE + G_SIZE + OFFSET_SIZE + ROTATION_SIZE;
    static const int ACCEL_BIAS_OFFSET = 0;
    static const int GYRO_BIAS_OFFSET = ACCEL_BIAS_SIZE;
    static const int G_OFFSET = 0;
    static const int TRANS_OFFSET = G_OFFSET + G_SIZE;
    static const int ROTATION_OFFSET = TRANS_OFFSET + OFFSET_SIZE;

    SparseHashAccumulator accum;
    double error;
    size_t bias_block_offset, gravity_block_offset, opt_size;

    SplineState* spline;
    CalibParam* calib_param;
    const Parameters* param;

    const bool pose_fixed;

    Linearizer(size_t _bias_block_offset, size_t _gravity_block_offset, size_t _opt_size, SplineState* spl,
        CalibParam* cbpar, const Parameters* par, const bool _pose_fixed)
        : bias_block_offset(_bias_block_offset), gravity_block_offset(_gravity_block_offset), opt_size(_opt_size),
          spline(spl), calib_param(cbpar), param(par),pose_fixed(_pose_fixed)
    {
        accum.reset(opt_size);
        error = 0;
    }

    ~Linearizer() {}

    void operator()(const Eigen::aligned_deque<ImuData>& r) {
        size_t set_fixed = 1;
        Eigen::Vector3d accel_var_inv = param->accel_var_inv;
        Eigen::Vector3d gyro_var_inv = param->gyro_var_inv;
        const double w_acc = param->w_acc;
        const double w_gyro = param->w_gyro;
        accel_var_inv *= w_acc;
        gyro_var_inv *= w_gyro;
        accel_var_inv = accel_var_inv.cwiseProduct(accel_var_inv);
        gyro_var_inv = gyro_var_inv.cwiseProduct(gyro_var_inv);
        Eigen::Vector3d bias_accel_var_inv = param->bias_accel_var_inv;
        Eigen::Vector3d bias_gyro_var_inv = param->bias_gyro_var_inv;
        bias_accel_var_inv *= param->w_bias_acc;
        bias_gyro_var_inv *= param->w_bias_gyro;
        bias_accel_var_inv = bias_accel_var_inv.cwiseProduct(bias_accel_var_inv);
        bias_gyro_var_inv = bias_gyro_var_inv.cwiseProduct(bias_gyro_var_inv);
        double num_imu = r.size();
        accel_var_inv /= num_imu;
        gyro_var_inv /= num_imu;
        bias_accel_var_inv /= (num_imu - 1);
        bias_gyro_var_inv /= (num_imu - 1);
        for (const auto& pm : r) {
            Jacobian36 J_accel;
            Jacobian33 J_gyro;
            Jacobian J_bias;
            Eigen::Matrix3d J_bias_a, J_bias_g;
            Eigen::Matrix<double, 3, 2> J_g;
            int64_t t = pm.time_ns;
            Eigen::Matrix<double, 6, 1> residual;
            residual = Residuals::imuResidualJacobian(t, spline, &pm.accel, &pm.gyro, calib_param->gravity,
                                                      &J_accel, &J_gyro, &J_bias, &J_g);
            const Eigen::Vector3d r_a = residual.segment<3>(3);
            error += r_a.transpose() * accel_var_inv.asDiagonal() * r_a;
            size_t start_g = gravity_block_offset;
            size_t num_Ji = J_accel.d_val_d_knot.size();
            for (size_t i = 0; i < num_Ji; i++) {
                size_t start_i = (J_accel.start_idx + i) * POSE_SIZE;
                if (pose_fixed && start_i < set_fixed * POSE_SIZE) {
                    continue;
                }
                for (size_t j = 0; j <= i; j++) {
                    size_t start_j = (J_accel.start_idx + j) * POSE_SIZE;
                    if (pose_fixed && start_j < set_fixed * POSE_SIZE) {
                        continue;
                    }
                    accum.addH<POSE_SIZE, POSE_SIZE>(start_i, start_j,
                    J_accel.d_val_d_knot[i].transpose() * accel_var_inv.asDiagonal() * J_accel.d_val_d_knot[j]);
                }
                for (size_t j = 0; j < num_Ji; j++) {
                    size_t start_bias_a = bias_block_offset + (J_bias.start_idx + j) * BIAS_SIZE + ACCEL_BIAS_OFFSET;
                    accum.addH<ACCEL_BIAS_SIZE, POSE_SIZE>(start_bias_a, start_i,
                        J_bias.d_val_d_knot[j] * accel_var_inv.asDiagonal() * J_accel.d_val_d_knot[i]);
                }
                if (param->if_opt_g) {
                    accum.addH<G_SIZE, POSE_SIZE>(start_g, start_i,
                        J_g.transpose() * accel_var_inv.asDiagonal() * J_accel.d_val_d_knot[i]);
                }
                accum.addB<POSE_SIZE>(start_i, J_accel.d_val_d_knot[i].transpose() * accel_var_inv.asDiagonal() * r_a);
            }
            for (size_t i = 0; i < num_Ji; i++) {
                size_t start_bias_ai = bias_block_offset + (J_bias.start_idx + i) * BIAS_SIZE + ACCEL_BIAS_OFFSET;
                for (size_t j = 0; j <= i; j++) {
                    size_t start_bias_aj = bias_block_offset + (J_bias.start_idx + j) * BIAS_SIZE + ACCEL_BIAS_OFFSET;
                    Eigen::Matrix3d JT_w_J = J_bias.d_val_d_knot[i] * accel_var_inv.asDiagonal() * J_bias.d_val_d_knot[j];
                    accum.addH<ACCEL_BIAS_SIZE, ACCEL_BIAS_SIZE>(start_bias_ai, start_bias_aj, JT_w_J);
                }
                Eigen::Vector3d JT_w_r = J_bias.d_val_d_knot[i] * accel_var_inv.asDiagonal() * r_a;
                accum.addB<ACCEL_BIAS_SIZE>(start_bias_ai, JT_w_r);
                if (param->if_opt_g) {
                    accum.addH<G_SIZE, ACCEL_BIAS_SIZE>(start_g, start_bias_ai, J_g.transpose() *
                                                                 accel_var_inv.asDiagonal() * J_bias.d_val_d_knot[i]);
                }
            }
            if (param->if_opt_g) {
                accum.addH<G_SIZE, G_SIZE>(start_g, start_g, J_g.transpose() * accel_var_inv.asDiagonal() * J_g);
                accum.addB<G_SIZE>(start_g, J_g.transpose() * accel_var_inv.asDiagonal() * r_a);
            }
            const Eigen::Vector3d r_g = residual.head(3);
            error += r_g.transpose() * gyro_var_inv.asDiagonal() * r_g;
            for (size_t i = 0; i < num_Ji; i++) {
                size_t start_i = (J_gyro.start_idx + i) * POSE_SIZE + ROT_OFFSET;
                if (pose_fixed && start_i < set_fixed * POSE_SIZE) {
                    continue;
                }
                for (size_t j = 0; j <= i; j++) {
                    size_t start_j = (J_gyro.start_idx + j) * POSE_SIZE + ROT_OFFSET;
                    if (pose_fixed && start_j < set_fixed * POSE_SIZE) {
                        continue;
                    }
                accum.addH<ROT_SIZE, ROT_SIZE>(start_i, start_j, J_gyro.d_val_d_knot[i].transpose() *
                                                        gyro_var_inv.asDiagonal() * J_gyro.d_val_d_knot[j]);
                }
                for (size_t j = 0; j < num_Ji; j++) {
                    size_t start_bias_g = bias_block_offset + (J_bias.start_idx + j) * BIAS_SIZE + GYRO_BIAS_OFFSET;
                    accum.addH<GYRO_BIAS_SIZE, ROT_SIZE>(start_bias_g, start_i, J_bias.d_val_d_knot[j] *
                        gyro_var_inv.asDiagonal() * J_gyro.d_val_d_knot[i]);
                }
                accum.addB<ROT_SIZE>(start_i, J_gyro.d_val_d_knot[i].transpose() * gyro_var_inv.asDiagonal() * r_g);
            }
            for (size_t i = 0; i < num_Ji; i++) {
                size_t start_bias_gi = bias_block_offset + (J_bias.start_idx + i) * BIAS_SIZE + GYRO_BIAS_OFFSET;

                for (size_t j = 0; j <= i; j++) {
                    size_t start_bias_gj = bias_block_offset + (J_bias.start_idx + j) * BIAS_SIZE + GYRO_BIAS_OFFSET;

                    Eigen::Matrix3d JT_w_J = J_bias.d_val_d_knot[i] * gyro_var_inv.asDiagonal() * J_bias.d_val_d_knot[j];
                    accum.addH<GYRO_BIAS_SIZE, GYRO_BIAS_SIZE>(start_bias_gi, start_bias_gj, JT_w_J);
                }
                Eigen::Vector3d JT_w_r = J_bias.d_val_d_knot[i] * gyro_var_inv.asDiagonal() * r_g;
                accum.addB<GYRO_BIAS_SIZE>(start_bias_gi, JT_w_r);
            }
        }
        Eigen::aligned_deque<ImuData>::const_iterator it = r.begin();
        Jacobian Jb0;
        Eigen::Matrix<double, 6, 1> b0 = spline->itpBias((*it).time_ns, &Jb0);
        int64_t t0 = (*it).time_ns;
        it++;
        size_t num_J0 = Jb0.d_val_d_knot.size();
        while (it != r.end()) {
            Jacobian Jb1;
            Eigen::Matrix<double, 6, 1> b1 = spline->itpBias((*it).time_ns, &Jb1);
            int64_t t1 = (*it).time_ns;
            size_t num_J1 = Jb1.d_val_d_knot.size();
            Eigen::Vector3d r_ba = b1.head<3>() - b0.head<3>();
            Eigen::Vector3d r_bg = b1.tail<3>() - b0.tail<3>();
            error += r_ba.transpose() * bias_accel_var_inv.asDiagonal() * r_ba;
            error += r_bg.transpose() * bias_gyro_var_inv.asDiagonal() * r_bg;
            size_t delta_idx = Jb1.start_idx - Jb0.start_idx;
            delta_idx = delta_idx > 4 ? 4 : delta_idx;
            size_t max_num_cp = std::max(num_J0, num_J1);
            Eigen::aligned_vector<std::pair<size_t, double>> vJb(max_num_cp + delta_idx);
            for (size_t i = 0; i < max_num_cp; i++) {
                bool set_idx = false;
                if (i < num_J0) {
                    vJb[i].first = Jb0.start_idx + i;
                    set_idx = true;
                    vJb[i].second = - Jb0.d_val_d_knot[i];
                }
                if (i >= delta_idx) {
                    if (!set_idx)
                        vJb[i].first = Jb1.start_idx + i - delta_idx;
                    vJb[i].second += Jb1.d_val_d_knot[i - delta_idx];
                }
            }
            for (size_t i = 0; i < delta_idx; i++) {
                vJb[i + max_num_cp].first = Jb1.start_idx + i + max_num_cp - delta_idx;
                vJb[i + max_num_cp].second = Jb1.d_val_d_knot[max_num_cp - delta_idx + i];
            }
            for (size_t i = 0; i < vJb.size(); i++) {
                size_t start_i = bias_block_offset + vJb[i].first * BIAS_SIZE;
                size_t start_bias_ai = start_i  + ACCEL_BIAS_OFFSET;
                size_t start_bias_gi = start_i  + GYRO_BIAS_OFFSET;
                for (size_t j = 0; j <= i; j++) {
                    size_t start_j = bias_block_offset + vJb[j].first * BIAS_SIZE;
                    size_t start_bias_aj = start_j  + ACCEL_BIAS_OFFSET;
                    size_t start_bias_gj = start_j  + GYRO_BIAS_OFFSET;
                    double JT_J = vJb[i].second * vJb[j].second;
                    Eigen::Matrix3d JT_wba_J = JT_J * bias_accel_var_inv.asDiagonal();
                    Eigen::Matrix3d JT_wbg_J = JT_J * bias_gyro_var_inv.asDiagonal();
                    accum.addH<ACCEL_BIAS_SIZE, ACCEL_BIAS_SIZE>(start_bias_ai, start_bias_aj, JT_wba_J);
                    accum.addH<GYRO_BIAS_SIZE, GYRO_BIAS_SIZE>(start_bias_gi, start_bias_gj, JT_wbg_J);
                }
                accum.addB<ACCEL_BIAS_SIZE>(start_bias_ai, vJb[i].second * bias_accel_var_inv.asDiagonal() * r_ba);
                accum.addB<GYRO_BIAS_SIZE>(start_bias_gi, vJb[i].second * bias_gyro_var_inv.asDiagonal() * r_bg);
            }
            b0 = b1;
            Jb0 = Jb1;
            t0 = t1;
            num_J0 = num_J1;
            it++;
        }
    }

    void operator()(Eigen::aligned_deque<TOAData>& r)
    {
        size_t set_fixed = 1;
        bool if_reject_uwb_in_optimization = param->if_reject_uwb_in_optimization;
        double reject_uwb_thresh = param->reject_uwb_thresh;
        Eigen::aligned_map<uint16_t, double> time_offset = param->toa_offset;
        double w_uwb = param->w_uwb;
        double num_toa = r.size();
        Eigen::aligned_map<uint16_t, Eigen::Vector3d> an_list = param->anchor_list;
        Eigen::Vector3d t_UW = calib_param->t_nav_uwb;
        Eigen::Quaterniond q_UW = calib_param->q_nav_uwb;
        for (auto& pm : r) {
            Jacobian16 J;
            Eigen::Vector3d J_tUW, J_qUW;
            int64_t t_ns = pm.time_ns;
            uint16_t anchor_id = pm.anchor_id;
            Eigen::Vector3d offset = calib_param->offset;
            double residual = Residuals::toaResidualJacobian(t_ns, spline, pm.data, an_list.at(anchor_id), offset,
                                                             t_UW, q_UW, time_offset[anchor_id], &J, &J_tUW, &J_qUW);
            if (if_reject_uwb_in_optimization && abs(residual) > reject_uwb_thresh) {
                continue;
            }
            double e2 = residual * residual;
            double range_std_inv = 1.0 / sqrt(num_toa);
            size_t num_J = J.d_val_d_knot.size();
            error += e2 * range_std_inv * range_std_inv * w_uwb * w_uwb;
            for (size_t i = 0; i < num_J; i++) {
                J.d_val_d_knot[i] *= range_std_inv;
            }
            J_tUW *= range_std_inv;
            J_qUW *= range_std_inv;
            residual *= range_std_inv;
            for (size_t i = 0; i < num_J; i++) {
                J.d_val_d_knot[i] *= w_uwb;
            }
            J_tUW *= w_uwb;
            J_qUW *= w_uwb;
            residual *= w_uwb;
            size_t start_tUW = gravity_block_offset + TRANS_OFFSET;
            size_t start_qUW = gravity_block_offset + ROTATION_OFFSET;
            for (size_t i = 0; i < num_J; i++) {
                size_t start_i = (J.start_idx + i) * POSE_SIZE;
                if (pose_fixed && start_i < set_fixed * POSE_SIZE) {
                    continue;
                }
                for (size_t j = 0; j <= i; j++) {
                    size_t start_j = (J.start_idx + j) * POSE_SIZE;
                    if (pose_fixed && start_j < set_fixed * POSE_SIZE) {
                        continue;
                    }
                    accum.addH<POSE_SIZE, POSE_SIZE>(start_i, start_j, J.d_val_d_knot[i] * J.d_val_d_knot[j].transpose());
                }
                if (param->if_opt_transform) {
                    accum.addH<OFFSET_SIZE, POSE_SIZE>(start_tUW, start_i, J_tUW * J.d_val_d_knot[i].transpose());
                    accum.addH<ROTATION_SIZE, POSE_SIZE>(start_qUW, start_i, J_qUW * J.d_val_d_knot[i].transpose());
                }
                accum.addB<POSE_SIZE>(start_i, J.d_val_d_knot[i] * residual);
            }
            if (param->if_opt_transform) {
                accum.addH<OFFSET_SIZE, OFFSET_SIZE>(start_tUW, start_tUW, J_tUW * J_tUW.transpose());
                accum.addB<OFFSET_SIZE>(start_tUW, J_tUW * residual);
                accum.addH<ROTATION_SIZE, OFFSET_SIZE>(start_qUW, start_tUW, J_qUW * J_tUW.transpose());
                accum.addH<ROTATION_SIZE, ROTATION_SIZE>(start_qUW, start_qUW, J_qUW * J_qUW.transpose());
                accum.addB<ROTATION_SIZE>(start_qUW, J_qUW * residual);
            }
        }
    }

    void operator()(Eigen::aligned_deque<TDOAData>& r)
    {
        size_t set_fixed = 1;
        bool if_reject_uwb_in_optimization = param->if_reject_uwb_in_optimization;
        double reject_uwb_thresh = param->reject_uwb_thresh;
        double w_uwb = param->w_uwb;
        double num_tdoa = r.size();
        Eigen::aligned_map<uint16_t, Eigen::Vector3d> an_list = param->anchor_list;
        Eigen::Vector3d t_UW = calib_param->t_nav_uwb;
        Eigen::Quaterniond q_UW = calib_param->q_nav_uwb;
        for (auto& pm : r) {
            Jacobian16 J;
            Eigen::Vector3d J_tUW, J_qUW;
            int64_t t_ns = pm.time_ns;
            int anchor_idA = pm.idA;
            int anchor_idB = pm.idB;
            Eigen::Vector3d offset = calib_param->offset;
            double residual = Residuals::tdoaResidualJacobian(t_ns, spline, pm.data, an_list.at(anchor_idA),
                                                              an_list.at(anchor_idB), offset, t_UW, q_UW, &J, &J_tUW, &J_qUW);
            size_t num_Ji = J.d_val_d_knot.size();
            if (if_reject_uwb_in_optimization && abs(residual) > reject_uwb_thresh) {
              continue;
            }
            double range_std_inv = 1.0 / sqrt(num_tdoa);
            for (size_t i = 0; i < num_Ji; i++) {
                J.d_val_d_knot[i] *= range_std_inv;
            }
            J_tUW *= range_std_inv;
            J_qUW *= range_std_inv;
            residual *= range_std_inv;
            for (size_t i = 0; i < num_Ji; i++) {
                J.d_val_d_knot[i] *= w_uwb;
            }
            J_tUW *= w_uwb;
            J_qUW *= w_uwb;
            residual *= w_uwb;
            error += residual * residual;
            size_t start_tUW = gravity_block_offset + TRANS_OFFSET;
            size_t start_qUW = gravity_block_offset + ROTATION_OFFSET;
            for (size_t i = 0; i < num_Ji; i++) {
                size_t start_i = (J.start_idx + i) * POSE_SIZE;
                if (pose_fixed && start_i < set_fixed * POSE_SIZE) {
                    continue;
                }
                for (size_t j = 0; j <= i; j++) {
                    size_t start_j = (J.start_idx + j) * POSE_SIZE;
                    if (pose_fixed && start_j < set_fixed * POSE_SIZE) {
                        continue;
                    }
                    accum.addH<POSE_SIZE, POSE_SIZE>(start_i, start_j, J.d_val_d_knot[i] * J.d_val_d_knot[j].transpose());
                }
                if (param->if_opt_transform) {
                    accum.addH<OFFSET_SIZE, POSE_SIZE>(start_tUW, start_i, J_tUW * J.d_val_d_knot[i].transpose());
                    accum.addH<ROTATION_SIZE, POSE_SIZE>(start_qUW, start_i, J_qUW * J.d_val_d_knot[i].transpose());
                }
                accum.addB<POSE_SIZE>(start_i, J.d_val_d_knot[i] * residual);
            }
            if (param->if_opt_transform) {
                accum.addH<OFFSET_SIZE, OFFSET_SIZE>(start_tUW, start_tUW, J_tUW * J_tUW.transpose());
                accum.addB<OFFSET_SIZE>(start_tUW, J_tUW * residual);
                accum.addH<ROTATION_SIZE, OFFSET_SIZE>(start_qUW, start_tUW, J_qUW * J_tUW.transpose());
                accum.addH<ROTATION_SIZE, ROTATION_SIZE>(start_qUW, start_qUW, J_qUW * J_qUW.transpose());
                accum.addB<ROTATION_SIZE>(start_qUW, J_qUW * residual);
            }
        }
     }
     EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct ComputeErrorSplineOpt
{
    double error;
    SplineState* spline;
    CalibParam* calib_param;
    const Parameters* param;

    ComputeErrorSplineOpt(SplineState* spl, CalibParam* cbpar, const Parameters* par)
        : spline(spl), calib_param(cbpar), param(par)
    {
        error = 0;
    }

    ~ComputeErrorSplineOpt() {}

    void operator()(const Eigen::aligned_deque<ImuData>& r)
    {
        Eigen::Vector3d accel_var_inv = param->accel_var_inv;
        Eigen::Vector3d gyro_var_inv = param->gyro_var_inv;
        const double w_acc = param->w_acc;
        const double w_gyro = param->w_gyro;
        accel_var_inv *= w_acc;
        gyro_var_inv *= w_gyro;
        accel_var_inv = accel_var_inv.cwiseProduct(accel_var_inv);
        gyro_var_inv = gyro_var_inv.cwiseProduct(gyro_var_inv);
        double num_imu = r.size();
        accel_var_inv /= num_imu;
        gyro_var_inv /= num_imu;
        Eigen::Vector3d bias_accel_var_inv = param->bias_accel_var_inv;
        Eigen::Vector3d bias_gyro_var_inv = param->bias_gyro_var_inv;
        bias_accel_var_inv *= param->w_bias_acc;
        bias_gyro_var_inv *= param->w_bias_gyro;
        bias_accel_var_inv = bias_accel_var_inv.cwiseProduct(bias_accel_var_inv);
        bias_gyro_var_inv = bias_gyro_var_inv.cwiseProduct(bias_gyro_var_inv);
        bias_accel_var_inv /= (num_imu - 1);
        bias_gyro_var_inv /= (num_imu - 1);
        for (const auto& pm : r) {
            int64_t t = pm.time_ns;
            Eigen::Matrix<double, 6, 1> residual = Residuals::imuResidual(t, spline, &pm.accel, &pm.gyro, calib_param->gravity);
            const Eigen::Vector3d r_a = residual.segment<3>(3);
            error += r_a.transpose() * accel_var_inv.asDiagonal() * r_a;
            const Eigen::Vector3d r_g = residual.head(3);
            error += r_g.transpose() * gyro_var_inv.asDiagonal() * r_g;
        }
        Eigen::aligned_deque<ImuData>::const_iterator it = r.begin();
        Eigen::Matrix<double, 6, 1> b0 = spline->itpBias((*it).time_ns);
        it ++;
        while (it != r.end()) {
            Eigen::Matrix<double, 6, 1> b1 = spline->itpBias((*it).time_ns);
            Eigen::Vector3d r_ba = b1.head<3>() - b0.head<3>();
            Eigen::Vector3d r_bg = b1.tail<3>() - b0.tail<3>();
            error += r_ba.transpose() * bias_accel_var_inv.asDiagonal() * r_ba;
            error += r_bg.transpose() * bias_gyro_var_inv.asDiagonal() * r_bg;
            b0 = b1;
            it++;
        }
    }

    void operator() (const Eigen::aligned_deque<TOAData>& r)
    {
        bool if_reject_uwb_in_optimization = param->if_reject_uwb_in_optimization;
        double reject_uwb_thresh = param->reject_uwb_thresh;
        Eigen::aligned_map<uint16_t, double> time_offset = param->toa_offset;
        double w_uwb = param->w_uwb;
        double num_toa = r.size();
        Eigen::aligned_map<uint16_t, Eigen::Vector3d> an_list = param->anchor_list;
        Eigen::Vector3d t_UW = calib_param->t_nav_uwb;
        Eigen::Quaterniond q_UW = calib_param->q_nav_uwb;
        for (const auto& pm : r) {
            int64_t t_ns = pm.time_ns;
            Eigen::Vector3d offset = calib_param->offset;
            uint16_t anchor_id = pm.anchor_id;
            double residual = Residuals::toaResidual(t_ns, spline, pm.data, an_list.at(anchor_id),
                                                     offset, t_UW, q_UW, time_offset[anchor_id]);
            if (if_reject_uwb_in_optimization && abs(residual) > reject_uwb_thresh) {
              continue;
            }
            double e2 = residual * residual;
            double range_std_inv = 1.0 / sqrt(num_toa);
            error += e2 * range_std_inv * range_std_inv * w_uwb * w_uwb;
        }
    }

    void operator()(Eigen::aligned_deque<TDOAData>& r)
    {
        bool if_reject_uwb_in_optimization = param->if_reject_uwb_in_optimization;
        double reject_uwb_thresh = param->reject_uwb_thresh;
        double w_uwb = param->w_uwb;
        double num_tdoa = r.size();
        Eigen::aligned_map<uint16_t, Eigen::Vector3d> an_list = param->anchor_list;
        Eigen::Vector3d t_UW = calib_param->t_nav_uwb;
        Eigen::Quaterniond q_UW = calib_param->q_nav_uwb;
        for (auto& pm : r) {
            Jacobian16 J;
            Eigen::Vector3d J_tUW, J_qUW;
            int64_t t_ns = pm.time_ns;
            Eigen::Vector3d offset = calib_param->offset;
            double residual = Residuals::tdoaResidual(t_ns, spline, pm.data, an_list.at(pm.idA), an_list.at(pm.idB),
                                                      offset, t_UW, q_UW);
            if (if_reject_uwb_in_optimization && abs(residual) > reject_uwb_thresh) {
                continue;
            }
            double range_std_inv = 1.0 / sqrt(num_tdoa);
            double e2 = residual * residual;
            error += e2 * w_uwb * w_uwb * range_std_inv * range_std_inv;
        }
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
