#pragma once

#include "utils/math_tools.h"
#include "SplineState.h"

class Residuals
{

public:

    static Eigen::Matrix<double, 6, 1> imuResidual(int64_t time_ns, const SplineState* spl,
       const Eigen::Vector3d* acc,  const Eigen::Vector3d* gyro, const Eigen::Vector3d& g)
    {
        Eigen::Quaterniond q_itp;
        Eigen::Vector3d rot_vel;
        spl->itpQuaternion(time_ns, &q_itp, &rot_vel);
        Eigen::Matrix<double,6, 1> bias = spl->itpBias(time_ns);
        Eigen::Vector3d a_w = spl->itpPosition<2>(time_ns) + g;
        Eigen::Matrix3d RT = q_itp.inverse().toRotationMatrix();
        Eigen::Vector3d acc_res = RT * a_w - *acc + bias.head<3>();
        Eigen::Vector3d gyro_res = rot_vel - *gyro + bias.tail<3>();
        Eigen::Matrix<double, 6, 1> residual;
        residual.head<3>() = gyro_res;
        residual.tail<3>() = acc_res;
        return residual;
    }

    static Eigen::Matrix<double, 6, 1> imuResidualJacobian(int64_t time_ns, const SplineState* spl,
        const Eigen::Vector3d* acc,  const Eigen::Vector3d* gyro, const Eigen::Vector3d& g,
        Jacobian36* J_acc, Jacobian33* J_gyro, Jacobian* J_bias = nullptr, Eigen::Matrix<double, 3, 2>* J_g = nullptr)
    {
        Eigen::Quaterniond q_itp;
        Eigen::Vector3d rot_vel;
        Jacobian43 J_rot;
        Jacobian J_line_acc;
        spl->itpQuaternion(time_ns, &q_itp, &rot_vel, &J_rot, J_gyro);
        Eigen::Vector3d a_w_no_g = spl->itpPosition<2>(time_ns, &J_line_acc);
        Eigen::Matrix<double,6, 1> bias = spl->itpBias(time_ns, J_bias);
        Eigen::Vector3d a_w = a_w_no_g + g;
        Eigen::Matrix3d RT = q_itp.inverse().toRotationMatrix();
        Eigen::Vector3d acc_res = RT * a_w - *acc + bias.head<3>();
        Eigen::Vector3d gyro_res = rot_vel - *gyro + bias.tail<3>();
        Eigen::Matrix<double, 6, 1> residual;
        residual.head<3>() = gyro_res;
        residual.tail<3>() = acc_res;
        Eigen::Matrix<double, 3, 4> tmp;
        Quater::drot(a_w, q_itp, tmp);
        J_acc->start_idx = J_line_acc.start_idx;
        int num_J = J_line_acc.d_val_d_knot.size();
        J_acc->d_val_d_knot.resize(num_J);
        for (int i = 0; i < num_J; i++) {
            J_acc->d_val_d_knot[i].topLeftCorner<3, 3>() = RT * J_line_acc.d_val_d_knot[i];
            J_acc->d_val_d_knot[i].bottomRightCorner<3, 3>() = tmp * J_rot.d_val_d_knot[i];
        }
        if (J_g) {
            (*J_g) = RT * Sphere::TangentBasis(g);
        }
        return residual;
    }

    static double tdoaResidual(int64_t time_ns, const SplineState* spl, const double tdoa,
        const Eigen::Vector3d& pos_anchor_i, const Eigen::Vector3d& pos_anchor_j,
        const Eigen::Vector3d& offset, const Eigen::Vector3d& t_UW, const Eigen::Quaterniond& q_UW)
    {
        Eigen::Quaterniond q_itp;
        spl->itpQuaternion(time_ns, &q_itp);
        Eigen::Matrix3d R_UW(q_UW);
        Eigen::Vector3d p_tag_W = q_itp * offset + spl->itpPosition(time_ns);
        Eigen::Vector3d p_tag_U = R_UW * p_tag_W + t_UW;
        Eigen::Vector3d diff_i = p_tag_U - pos_anchor_i;
        Eigen::Vector3d diff_j = p_tag_U - pos_anchor_j;
        return diff_j.norm() - diff_i.norm() - tdoa;
    }

    static double tdoaResidualJacobian(int64_t time_ns, const SplineState* spl, const double tdoa,
        const Eigen::Vector3d& pos_anchor_i, const Eigen::Vector3d& pos_anchor_j,
        const Eigen::Vector3d& offset, const Eigen::Vector3d& t_UW, const Eigen::Quaterniond& q_UW,
        Jacobian16 *J, Eigen::Vector3d *J_tUW = nullptr, Eigen::Vector3d *J_qUW = nullptr)
    {
        Jacobian J_pos;
        Jacobian43 J_quat;
        Eigen::Quaterniond q_itp;
        spl->itpQuaternion(time_ns, &q_itp, nullptr, &J_quat);
        Eigen::Matrix3d R_UW(q_UW);
        Eigen::Vector3d p_tag_W = q_itp * offset + spl->itpPosition(time_ns, &J_pos);
        Eigen::Vector3d p_tag_U = R_UW * p_tag_W + t_UW;
        Eigen::Vector3d diff_i = p_tag_U - pos_anchor_i;
        Eigen::Vector3d diff_j = p_tag_U - pos_anchor_j;
        Eigen::Vector3d dtdoa = diff_j.normalized() - diff_i.normalized();
        Eigen::Matrix<double, 3, 4> drot1, drot2;
        Quater::drot2(p_tag_W, q_UW, drot1);
        Quater::drot2(offset, q_itp, drot2);
        J->start_idx = J_pos.start_idx;
        int num_J = J_quat.d_val_d_knot.size();
        J->d_val_d_knot.resize(num_J);
        Eigen::Matrix<double, 1, 3> tmp = dtdoa.transpose() * R_UW;
        for (int i = 0; i < num_J; i++) {
            J->d_val_d_knot[i].head<3>() = tmp * J_pos.d_val_d_knot[i];
            J->d_val_d_knot[i].tail<3>() = tmp * drot2 * J_quat.d_val_d_knot[i];
        }
        if (J_tUW) {
            *J_tUW = dtdoa;
        }
        if (J_qUW) {         
            Eigen::Matrix4d delta_Q_l;
            Quater::Qleft(q_UW, delta_Q_l);
            *J_qUW = dtdoa.transpose() * drot1 * delta_Q_l.rightCols(3);
        }
        return diff_j.norm() - diff_i.norm() - tdoa;
    }

    static double toaResidual(int64_t time_ns, const SplineState* spl, const double toa,
        const Eigen::Vector3d& pos_anchor, const Eigen::Vector3d& offset, const Eigen::Vector3d& t_UW,
        const Eigen::Quaterniond& q_UW, const double distance_offset)
    {
        Eigen::Quaterniond q_itp;
        spl->itpQuaternion(time_ns, &q_itp);
        Eigen::Matrix3d R_UW(q_UW);
        Eigen::Vector3d W_p_tag = q_itp * offset + spl->itpPosition(time_ns);
        Eigen::Vector3d diff = R_UW * W_p_tag + t_UW - pos_anchor;
        return diff.norm() - toa - distance_offset;
    }

    static double toaResidualJacobian(int64_t time_ns, const SplineState* spl, const double toa,
        const Eigen::Vector3d& pos_anchor, const Eigen::Vector3d& offset, const Eigen::Vector3d& t_UW,
        const Eigen::Quaterniond& q_UW, const double distance_offset,
        Jacobian16 *J, Eigen::Vector3d *J_tUW = nullptr, Eigen::Vector3d *J_qUW = nullptr)
    {
        Jacobian J_pos;
        Jacobian43 J_quat;
        Eigen::Quaterniond q_itp;
        Eigen::Matrix<double, 3, 4> drot1, drot2;
        spl->itpQuaternion(time_ns, &q_itp, nullptr, &J_quat);
        Eigen::Matrix3d R_UW(q_UW);
        Eigen::Vector3d W_p_tag = q_itp * offset + spl->itpPosition(time_ns, &J_pos);
        Eigen::Vector3d diff = R_UW * W_p_tag + t_UW - pos_anchor;
        Quater::drot2(W_p_tag, q_UW, drot1);
        Quater::drot2((offset), q_itp, drot2);
        Eigen::Vector3d dtoa = diff.normalized();
        J->start_idx = J_pos.start_idx;
        int num_J = J_quat.d_val_d_knot.size();
        J->d_val_d_knot.resize(num_J);
        Eigen::Matrix<double, 1, 3> tmp2 = dtoa.transpose() * R_UW;
        for (int i = 0; i < num_J; i++) {
            J->d_val_d_knot[i].head<3>() = tmp2 * J_pos.d_val_d_knot[i];
            J->d_val_d_knot[i].tail<3>() = tmp2 * drot2 * J_quat.d_val_d_knot[i];
        }
        if (J_tUW) {
            *J_tUW = dtoa;
        }
        if (J_qUW) {
            Eigen::Matrix4d delta_Q_l;
            Quater::Qleft(q_UW, delta_Q_l);
            *J_qUW = dtoa.transpose() * drot1 * delta_Q_l.rightCols(3);
        }
        return diff.norm() - toa - distance_offset;
    }
};
