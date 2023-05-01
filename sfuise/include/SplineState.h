#pragma once

#include "utils/common_utils.h"
#include "utils/math_tools.h"
#include "sfuise_msgs/Spline.h"
#include "sfuise_msgs/Knot.h"

template <class MatT>
struct JacobianStruct {
    size_t start_idx;
    std::vector<MatT> d_val_d_knot;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

typedef JacobianStruct<double> Jacobian;
typedef JacobianStruct<Eigen::Matrix<double, 4, 3>> Jacobian43;
typedef JacobianStruct<Eigen::Matrix3d> Jacobian33;
typedef JacobianStruct<Eigen::Vector3d> Jacobian13;
typedef JacobianStruct<Eigen::Matrix<double, 6, 1>> Jacobian16;
typedef JacobianStruct<Eigen::Matrix<double, 3, 6>> Jacobian36;

class SplineState
{

  public:

    SplineState() {};

    void init(int64_t dt_ns_, int num_knot_, int64_t start_t_ns_, int start_i_ = 0)
    {
        if_first = true;
        dt_ns = dt_ns_;
        start_t_ns = start_t_ns_;
        num_knot = num_knot_;
        inv_dt = 1e9 / dt_ns;
        start_i = start_i_;
        pow_inv_dt[0] = 1.0;
        pow_inv_dt[1] = inv_dt;
        pow_inv_dt[2] = inv_dt * inv_dt;
        pow_inv_dt[3] = pow_inv_dt[2] * inv_dt;
        Eigen::Quaterniond q0 = Eigen::Quaterniond::Identity();
        Eigen::Vector3d t0 = Eigen::Vector3d::Zero();
        Eigen::Matrix<double, 6, 1> b0 = Eigen::Matrix<double, 6, 1>::Zero();
        q_idle = {q0, q0, q0};
        t_idle = {t0, t0, t0};
        b_idle = {b0, b0, b0};
    }

    void setOneStateKnot(int i, Eigen::Quaterniond q, Eigen::Vector3d pos, Eigen::Matrix<double, 6, 1> bias)
    {
        t_knots[i] = pos;
        q_knots[i] = q;
        b_knots[i] = bias;
    }

    void updateKnots(SplineState* other)
    {
        size_t num_knots = other->numKnots();
        for (size_t i = 0; i < num_knots; i++) {
            if (i + other->start_i < num_knot) {
                setOneStateKnot(i + other->start_i, other->q_knots[i],other->t_knots[i],other->b_knots[i]);
            } else {
                addOneStateKnot(other->q_knots[i],other->t_knots[i],other->b_knots[i]);
            }
        }
        q_idle = other->q_idle;
        t_idle = other->t_idle;
        b_idle = other->b_idle;
    }

    Eigen::Vector3d getIdlePos(int idx)
    {
        return t_idle[idx];
    }

    void setIdles(int idx, Eigen::Vector3d t, Eigen::Quaterniond q, Eigen::Matrix<double, 6, 1> b)
    {
         t_idle[idx] = t;
         q_idle[idx] = q;
         b_idle[idx] = b;
    }

    void updateBiasIdleFirstWindow()
    {
        if (!if_first) {
            return;
        } else {
            b_idle = {b_knots[0], b_knots[0], b_knots[0]};
        }
    }

    void addOneStateKnot(Eigen::Quaterniond q, Eigen::Vector3d pos, Eigen::Matrix<double, 6, 1> bias)
    {
        if (num_knot > 1) {
            Eigen::Quaterniond q0 = q_knots[num_knot - 1];
            Eigen::Quaterniond q1 = q;
            double dot_product = q0.dot(q1);
            if (dot_product < 0) {
               q = Eigen::Quaterniond(-q.w(), -q.x(), -q.y(), -q.z());
            }
        }
        if (abs(q.norm() - 1) > 1e-5) {
            q.normalize();
        }
        q_knots.push_back(q);
        t_knots.push_back(pos);
        b_knots.push_back(bias);
        num_knot++;
    }

    void checkQuaternionControlPoints()
    {
        if (num_knot > 1) {
            for (size_t i = 1; i < num_knot; i++) {
                Eigen::Quaterniond q1 = q_knots[i];
                if (q_knots[i - 1].dot(q1) < 0) {
                    q_knots[i] = Eigen::Quaterniond(-q1.w(), -q1.x(), -q1.y(), -q1.z());
                }
            }
        }
    }

    void removeOneOldState()
    {
        q_idle = {q_idle[1], q_idle[2], q_knots.front()};
        t_idle = {t_idle[1], t_idle[2], t_knots.front()};
        b_idle = {b_idle[1], b_idle[2], b_knots.front()};
        q_knots.pop_front();
        t_knots.pop_front();
        b_knots.pop_front();
        num_knot--;
        start_i++;
        start_t_ns += dt_ns;
        if_first = false;
    }

    void getAllStateKnots(Eigen::aligned_deque<Eigen::Vector3d>& knots_trans,
        Eigen::aligned_deque<Eigen::Quaterniond>& knots_rot,
        Eigen::aligned_deque<Eigen::Matrix<double,6, 1>>& knots_bias)
    {
        knots_trans = t_knots;
        knots_rot = q_knots;
        knots_bias = b_knots;
    }

    void setAllKnots(Eigen::aligned_deque<Eigen::Vector3d>& knots_trans,
        Eigen::aligned_deque<Eigen::Quaterniond>& knots_rot,
        Eigen::aligned_deque<Eigen::Matrix<double,6, 1>>& knots_bias)
    {
        t_knots = knots_trans;
        q_knots = knots_rot;
        b_knots = knots_bias;
        updateBiasIdleFirstWindow();
    }

    int64_t getKnotTimeNs(size_t i) const
    {
        return start_t_ns + i * dt_ns;
    }

    Eigen::Quaterniond getKnotOrt(size_t i)
    {
        return q_knots[i];
    }

    Eigen::Vector3d getKnotPos(size_t i)
    {
        return t_knots[i];
    }

    Eigen::Matrix<double, 6, 1> getKnotBias(size_t i)
    {
        return b_knots[i];
    }

    Eigen::Vector3d extrapolatePosKnot(size_t idx)
    {
        Eigen::Quaterniond last_ort = q_knots[num_knot - idx - 1];
        Eigen::Quaterniond cur_ort = q_knots[num_knot - idx];
        Eigen::Vector3d last_trans = t_knots[num_knot - idx - 1];
        Eigen::Vector3d cur_trans = t_knots[num_knot - idx];
        Eigen::Vector3d rel_trans = last_ort.inverse() * (cur_trans - last_trans);
        return cur_trans + cur_ort * rel_trans;
    }

    Eigen::Quaterniond extrapolateOrtKnot(size_t idx)
    {
        Eigen::Quaterniond last_ort = q_knots[num_knot - idx - 1];
        Eigen::Quaterniond cur_ort = q_knots[num_knot - idx];
        return cur_ort * last_ort.inverse() * cur_ort;
    }

    void applyPoseInc(int i, const Eigen::Matrix<double, 6, 1> &inc)
    {
        t_knots[i] += inc.head<3>();
        Eigen::Quaterniond q_inc;
        Quater::exp(inc.tail<3>(), q_inc);
        q_knots[i] *= q_inc;
    }

    void applyBiasInc(int i, const Eigen::Matrix<double, 6, 1>& inc)
    {
        b_knots[i] += inc;
    }

    int64_t maxTimeNs()
    {
        if (num_knot == 1) {
           return start_t_ns;
        }
        return start_t_ns + (num_knot - 1) * dt_ns - 1;
    }

    int64_t minTimeNs()
    {
        return start_t_ns + dt_ns * (!if_first ?  -1 : 0);
    }

    int64_t nextMaxTimeNs()
    {
        return start_t_ns + num_knot * dt_ns - 1;
    }

    size_t numKnots()
    {
        return num_knot;
    }

    template <int Derivative = 0>
    Eigen::Vector3d itpPosition (int64_t time_ns, Jacobian* J = nullptr) const
    {
        return itpEuclidean<Eigen::Vector3d, Derivative>(time_ns, t_idle, t_knots, J);
    }

    Eigen::Matrix<double, 6, 1> itpBias (int64_t time_ns, Jacobian* J = nullptr) const
    {
        return itpEuclidean<Eigen::Matrix<double, 6, 1>>(time_ns, b_idle, b_knots, J);
    }

    void itpQuaternion(int64_t t_ns, Eigen::Quaterniond* q_out = nullptr,
        Eigen::Vector3d* w_out = nullptr, Jacobian43* J_q = nullptr, Jacobian33* J_w = nullptr) const
    {
        double u;
        int64_t idx0;
        int idx_r;
        std::array<Eigen::Quaterniond, 4> cps;
        prepareInterpolation(t_ns, q_idle, q_knots, idx0, u, cps, idx_r);
        Eigen::Vector4d p;
        Eigen::Vector4d coeff;
        baseCoeffsWithTime<0>(p, u);
        coeff = cumulative_blending_matrix * p;
        Eigen::Quaterniond q_delta[3];
        q_delta[0] = cps[0].inverse() * cps[1];
        q_delta[1] = cps[1].inverse() * cps[2];
        q_delta[2] = cps[2].inverse() * cps[3];
        Eigen::Vector3d t_delta[3];
        Eigen::Vector3d t_delta_scale[3];
        Eigen::Quaterniond q_delta_scale[3];
        Eigen::Quaterniond q_itps[4];
        Eigen::Vector3d w_itps[4];
        Eigen::Vector4d dcoeff;
        if (J_q || J_w) {
            Eigen::Matrix<double, 3, 4> dlog_dq[3];
            Eigen::Matrix<double, 4, 3> dexp_dt[3];
            Quater::dlog(q_delta[0], t_delta[0], dlog_dq[0]);
            Quater::dlog(q_delta[1], t_delta[1], dlog_dq[1]);
            Quater::dlog(q_delta[2], t_delta[2], dlog_dq[2]);
            t_delta_scale[0] = t_delta[0] * coeff[1];
            t_delta_scale[1] = t_delta[1] * coeff[2];
            t_delta_scale[2] = t_delta[2] * coeff[3];
            Quater::dexp(t_delta_scale[0], q_delta_scale[0], dexp_dt[0]);
            Quater::dexp(t_delta_scale[1], q_delta_scale[1], dexp_dt[1]);
            Quater::dexp(t_delta_scale[2], q_delta_scale[2], dexp_dt[2]);
            int size_J = std::min(idx_r + 1, 4);
            Eigen::Matrix4d d_X_d_dj[3];
            Eigen::Matrix<double, 3, 4> d_r_d_dj[3];
            Eigen::Quaterniond q_r_all[4];
            q_r_all[3] = Eigen::Quaterniond::Identity();
            for (int i = 2; i >= 0; i-- ) {
                q_r_all[i] = q_delta_scale[i] * q_r_all[i+1];
            }
            Eigen::Matrix4d Q_l[size_J - 1];
            Eigen::Matrix4d Q_r[size_J - 1];
            for (int i = 0; i < size_J - 1; i++) {
                Quater::Qleft(q_delta[4 - size_J + i], Q_l[i]);
                Quater::Qright(q_delta[4 - size_J + i], Q_r[i]);
            }
            if (J_q) {
                q_itps[0] = cps[0];
                q_itps[1] = q_itps[0] * q_delta_scale[0];
                q_itps[2] = q_itps[1] * q_delta_scale[1];
                q_itps[3] = q_itps[2] * q_delta_scale[2];
                *q_out = q_itps[3];
                Eigen::Matrix4d Q_l_all[3];
                Quater::Qleft(q_itps[0], Q_l_all[0]);
                Quater::Qleft(q_itps[1], Q_l_all[1]);
                Quater::Qleft(q_itps[2], Q_l_all[2]);
                for (int i = 2; i >= 0; i--) {
                    Eigen::Matrix4d Q_r_all;
                    Quater::Qright(q_r_all[i+1], Q_r_all);
                    d_X_d_dj[i].noalias() = coeff[i + 1] * Q_r_all * Q_l_all[i] * dexp_dt[i] * dlog_dq[i];
                }
                J_q->d_val_d_knot.resize(size_J);
                for (int i = 0; i < size_J; i++) {
                    J_q->d_val_d_knot[i].setZero();
                }
                for (int i = 0; i < size_J - 1; i++) {
                    J_q->d_val_d_knot[i].noalias() -= d_X_d_dj[4 - size_J + i] * Q_r[i].rightCols(3);
                    J_q->d_val_d_knot[i + 1].noalias() += d_X_d_dj[4 - size_J + i] * Q_l[i].rightCols(3);
                }
                J_q->start_idx = idx0;
                if (size_J == 4) {
                    Eigen::Matrix4d Q_r_all;
                    Eigen::Matrix4d Q0_left;
                    Quater::Qright(q_r_all[0], Q_r_all);
                    Quater::Qleft(cps[0], Q0_left);
                    J_q->d_val_d_knot[0].noalias() += Q_r_all * Q0_left.rightCols(3);
                } else {
                    Eigen::Matrix4d Q_left;
                    Quater::Qleft(q_delta[3 - size_J], Q_left);
                    J_q->d_val_d_knot[0].noalias() += d_X_d_dj[3 - size_J] * Q_left.rightCols(3);
                }
            }
            if (J_w) {
                baseCoeffsWithTime<1>(p, u);
                dcoeff = inv_dt * cumulative_blending_matrix * p;
                w_itps[0].setZero();
                w_itps[1] = 2 * dcoeff[1] * t_delta[0];
                w_itps[2] = q_delta_scale[1].inverse() * w_itps[1] + 2 * dcoeff[2] * t_delta[1];
                w_itps[3] = q_delta_scale[2].inverse() * w_itps[2] + 2 * dcoeff[3] * t_delta[2];
                *w_out = w_itps[3];
                Eigen::Matrix<double, 3, 4> drot_dq[3];
                Quater::drot(w_itps[0], q_delta_scale[0], drot_dq[0]);
                Quater::drot(w_itps[1], q_delta_scale[1], drot_dq[1]);
                Quater::drot(w_itps[2], q_delta_scale[2], drot_dq[2]);
                for (int i = 2; i >= 0; i--) {
                    Eigen::Matrix3d d_vel_d_dj = coeff[i + 1] * drot_dq[i] * dexp_dt[i];
                    d_vel_d_dj.noalias() += 2 * dcoeff[i + 1] * Eigen::Matrix3d::Identity();
                    d_r_d_dj[i].noalias() = q_r_all[i+1].inverse().toRotationMatrix() * d_vel_d_dj * dlog_dq[i];
                }
                J_w->d_val_d_knot.resize(size_J);
                for (int i = 0; i < size_J; i++) {
                    J_w->d_val_d_knot[i].setZero();
                }
                for (int i = 0; i < size_J - 1; i++) {
                    J_w->d_val_d_knot[i].noalias() -= d_r_d_dj[4 - size_J + i] * Q_r[i].rightCols(3);
                    J_w->d_val_d_knot[i + 1].noalias() += d_r_d_dj[4 - size_J + i] * Q_l[i].rightCols(3);
                }
                J_w->start_idx = idx0;
                if (size_J != 4) {
                    Eigen::Matrix4d Q_left;
                    Quater::Qleft(q_delta[4 - size_J - 1], Q_left);
                    J_w->d_val_d_knot[0].noalias() += d_r_d_dj[3 - size_J] *  Q_left.rightCols(3);
                }
            }
        } else {
            Quater::log(q_delta[0], t_delta[0]);
            Quater::log(q_delta[1], t_delta[1]);
            Quater::log(q_delta[2], t_delta[2]);
            t_delta_scale[0] = t_delta[0] * coeff[1];
            t_delta_scale[1] = t_delta[1] * coeff[2];
            t_delta_scale[2] = t_delta[2] * coeff[3];
            Quater::exp(t_delta_scale[0], q_delta_scale[0]);
            Quater::exp(t_delta_scale[1], q_delta_scale[1]);
            Quater::exp(t_delta_scale[2], q_delta_scale[2]);
            if (q_out) {
                q_itps[0] = cps[0];
                q_itps[1] = q_itps[0] * q_delta_scale[0];
                q_itps[2] = q_itps[1] * q_delta_scale[1];
                q_itps[3] = q_itps[2] * q_delta_scale[2];
                q_itps[3].normalize();
                *q_out = q_itps[3];
            }
            if (w_out) {
                baseCoeffsWithTime<1>(p, u);
                dcoeff = inv_dt * cumulative_blending_matrix * p;

                w_itps[0].setZero();
                w_itps[1] = 2 * dcoeff[1] * t_delta[0];
                w_itps[2] = q_delta_scale[1].inverse() * w_itps[1] + 2 * dcoeff[2] * t_delta[1];
                w_itps[3] = q_delta_scale[2].inverse() * w_itps[2] + 2 * dcoeff[3] * t_delta[2];
                *w_out = w_itps[3];
            }
        }
    }

    void getSplineMsg(sfuise_msgs::Spline& spline_msg)
    {
        spline_msg.dt = dt_ns;
        spline_msg.start_t = start_t_ns;
        spline_msg.start_idx = start_i;
        for (size_t i = 0; i < num_knot; i++) {
            sfuise_msgs::Knot knot_msg;
            knot_msg.position.x = t_knots[i].x();
            knot_msg.position.y = t_knots[i].y();
            knot_msg.position.z = t_knots[i].z();
            knot_msg.orientation.w = q_knots[i].w();
            knot_msg.orientation.x = q_knots[i].x();
            knot_msg.orientation.y = q_knots[i].y();
            knot_msg.orientation.z = q_knots[i].z();
            Eigen::Matrix<double, 6, 1> bias = b_knots[i];
            knot_msg.bias_acc.x = bias[0];
            knot_msg.bias_acc.y = bias[1];
            knot_msg.bias_acc.z = bias[2];
            knot_msg.bias_gyro.x = bias[3];
            knot_msg.bias_gyro.y = bias[4];
            knot_msg.bias_gyro.z = bias[5];
            spline_msg.knots.push_back(knot_msg);
        }
        for (int i = 0; i < 3; i++) {
            sfuise_msgs::Knot idle_msg;
            idle_msg.position.x = t_idle[i].x();
            idle_msg.position.y = t_idle[i].y();
            idle_msg.position.z = t_idle[i].z();
            idle_msg.orientation.w = q_idle[i].w();
            idle_msg.orientation.x = q_idle[i].x();
            idle_msg.orientation.y = q_idle[i].y();
            idle_msg.orientation.z = q_idle[i].z();
            Eigen::Matrix<double, 6, 1> bias = b_idle[i];
            idle_msg.bias_acc.x = bias[0];
            idle_msg.bias_acc.y = bias[1];
            idle_msg.bias_acc.z = bias[2];
            idle_msg.bias_gyro.x = bias[3];
            idle_msg.bias_gyro.y = bias[4];
            idle_msg.bias_gyro.z = bias[5];
            spline_msg.idles.push_back(idle_msg);
        }
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  private:

    bool if_first;

    static constexpr double S_TO_NS = 1e9;
    static const Eigen::Matrix4d blending_matrix;
    static const Eigen::Matrix4d base_coefficients;
    static const Eigen::Matrix4d cumulative_blending_matrix;

    int64_t dt_ns;
    double inv_dt;
    std::array<double, 4> pow_inv_dt;
    int num_knot;
    int64_t start_i;
    int64_t start_t_ns;

    std::array<Eigen::Quaterniond, 3> q_idle;
    std::array<Eigen::Vector3d, 3> t_idle;
    std::array<Eigen::Matrix<double, 6, 1>, 3> b_idle;
    Eigen::aligned_deque<Eigen::Quaterniond> q_knots;
    Eigen::aligned_deque<Eigen::Vector3d> t_knots;
    Eigen::aligned_deque<Eigen::Matrix<double, 6, 1>> b_knots;

    template <typename _KnotT, int Derivative = 0>
    _KnotT itpEuclidean(int64_t t_ns, const std::array<_KnotT, 3>& knots_idle,
                        const Eigen::aligned_deque<_KnotT>& knots, Jacobian* J = nullptr) const
    {
        double u;
        int64_t idx0;
        int idx_r;
        std::array<_KnotT,4> cps;
        prepareInterpolation(t_ns, knots_idle, knots, idx0, u, cps, idx_r);
        Eigen::Vector4d p, coeff;
        baseCoeffsWithTime<Derivative>(p, u);
        coeff = pow_inv_dt[Derivative] * (blending_matrix * p);
        _KnotT res_out = coeff[0] * cps[0] + coeff[1] * cps[1] + coeff[2] * cps[2] + coeff[3] * cps[3];
        if (J) {
            int size_J = std::min(idx_r + 1, 4);
            J->d_val_d_knot.resize(size_J);
            for (int i = 0; i < size_J; i++) {
                J->d_val_d_knot[i] = coeff[4 - size_J + i];
            }
            J->start_idx = idx0;
        }
        return res_out;
    }

    template<typename _KnotT>
    void prepareInterpolation(int64_t t_ns, const std::array<_KnotT, 3>& knots_idle,
                              const Eigen::aligned_deque<_KnotT>& knots, int64_t& idx0, double& u,
                              std::array<_KnotT,4>& cps, int& idx_r) const
    {
        int64_t t_ns_rel = t_ns - start_t_ns;
        int idx_l = floor(double(t_ns_rel) / double(dt_ns));
        idx_r = idx_l + 1;
        idx0 = std::max(idx_l - 2, 0);
        for (int i = 0; i < 2 - idx_l; i++) {
            cps[i] = knots_idle[i + idx_l + 1];
        }
        int idx_window = std::max(0, 2 - idx_l);
        for (int i = 0; i < std::min(idx_l + 2, 4); i++) {
            cps[i + idx_window] = knots[idx0 + i];
        }
        u = (t_ns - start_t_ns - idx_l * dt_ns) / double(dt_ns);
    }

    template <int Derivative, class Derived>
    static void baseCoeffsWithTime(const Eigen::MatrixBase<Derived>& res_const, double t)
    {
        EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 4);
        Eigen::MatrixBase<Derived>& res = const_cast<Eigen::MatrixBase<Derived>&>(res_const);
        res.setZero();
        res[Derivative] = base_coefficients(Derivative, Derivative);
        double ti = t;
        for (int j = Derivative + 1; j < 4; j++) {
            res[j] = base_coefficients(Derivative, j) * ti;
            ti = ti * t;
        }
    }

    template <bool _Cumulative = false>
    static Eigen::Matrix4d computeBlendingMatrix()
    {
        Eigen::Matrix4d m;
        m.setZero();
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                double sum = 0;
                for (int s = j; s < 4; ++s) {
                    sum += std::pow(-1.0, s - j) * binomialCoefficient(4, s - j) *
                    std::pow(4 - s - 1.0, 4 - 1.0 - i);
                }
                m(j, i) = binomialCoefficient(3, 3 - i) * sum;
            }
        }
        if (_Cumulative) {
            for (int i = 0; i < 4; i++) {
                for (int j = i + 1; j < 4; j++) {
                    m.row(i) += m.row(j);
                }
            }
        }
        uint64_t factorial = 1;
        for (int i = 2; i < 4; ++i) {
            factorial *= i;
        }
        return m / factorial;
    }

    constexpr static inline uint64_t binomialCoefficient(uint64_t n, uint64_t k)
    {
        if (k > n) return 0;
        uint64_t r = 1;
        for (uint64_t d = 1; d <= k; ++d) {
            r *= n--;
            r /= d;
        }
        return r;
    }

    static Eigen::Matrix4d computeBaseCoefficients()
    {
        Eigen::Matrix4d base_coeff;
        base_coeff.setZero();
        base_coeff.row(0).setOnes();
        int order = 3;
        for (int n = 1; n < 4; n++) {
            for (int i = 3 - order; i < 4; i++) {
                base_coeff(n, i) = (order - 3 + i) * base_coeff(n - 1, i);
            }
            order--;
        }
        return base_coeff;
    }
};

const Eigen::Matrix4d SplineState::base_coefficients = SplineState::computeBaseCoefficients();
const Eigen::Matrix4d SplineState::blending_matrix = SplineState::computeBlendingMatrix();
const Eigen::Matrix4d SplineState::cumulative_blending_matrix = SplineState::computeBlendingMatrix<true>();
