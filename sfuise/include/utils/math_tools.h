#pragma once

#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>
#include <boost/math/special_functions/sinc.hpp>
#include <cmath>

class Quater
{
  public:
    static void Qleft(const Eigen::Quaterniond &q, Eigen::Matrix4d &ans)
    {
        double q1 = q.w();
        double q2 = q.x();
        double q3 = q.y();
        double q4 = q.z();
        ans << q1, -q2, -q3, -q4,
               q2,  q1, -q4,  q3,
               q3,  q4,  q1, -q2,
               q4, -q3,  q2,  q1;
    }

    static void Qright(const Eigen::Quaterniond &q, Eigen::Matrix4d &ans)
    {
        double q1 = q.w();
        double q2 = q.x();
        double q3 = q.y();
        double q4 = q.z();
        ans << q1, -q2, -q3, -q4,
               q2,  q1,  q4, -q3,
               q3, -q4,  q1,  q2,
               q4,  q3, -q2,  q1;
    }

    static void exp(const Eigen::Vector3d &v, Eigen::Quaterniond& q)
    {
        double v_norm = v.norm();
        q.w() = std::cos(v_norm);
        q.vec() = boost::math::sinc_pi(v_norm) * v;
    }

    static void log(const Eigen::Quaterniond &q, Eigen::Vector3d& v)
    {
        Eigen::Quaterniond qn = q.normalized();
        Eigen::Vector3d rv = qn.vec();
        double norm_rv = rv.norm();
        if (norm_rv > 1e-5) {
            v = std::atan(norm_rv / qn.w()) * rv / norm_rv;
        } else {
            v.setZero();
        }
    }


    static void dexp(const Eigen::Vector3d& v, Eigen::Quaterniond& q, Eigen::Matrix<double,4,3>& J)
    {
        double v_norm = v.norm();
        if (v_norm == 0) {
            J.row(0).setZero();
            J.bottomRows<3>().setIdentity();
            q.setIdentity();
            return;
        }
        double sinc = boost::math::sinc_pi(v_norm);
        q.w() = std::cos(v_norm);
        q.vec() = sinc * v;
        double v1 = v(0);
        double v2 = v(1);
        double v3 = v(2);
        double tmp = (q.w() - sinc) / (v_norm * v_norm);
        double tmp_v1 = tmp * v1;
        double tmp_v2 = tmp * v2;
        double tmp_v11 = tmp_v1 * v1 + sinc;
        double tmp_v12 = tmp_v1 * v2;
        double tmp_v13 = tmp_v1 * v3;
        double tmp_v22 = tmp_v2 * v2 + sinc;
        double tmp_v23 = tmp_v2 * v3;
        double tmp_v33 = tmp * v3 * v3 + sinc;
        J << -v1 * sinc, -v2 * sinc, -v3 * sinc,
             tmp_v11, tmp_v12, tmp_v13,
             tmp_v12, tmp_v22, tmp_v23,
             tmp_v13, tmp_v23, tmp_v33;

    }

    static void dlog(const Eigen::Quaterniond &r, Eigen::Vector3d& v, Eigen::Matrix<double, 3, 4>& J)
    {
        Eigen::Quaterniond rn = r.normalized();
        double rnw = rn.w();
        double r1 = rn.x();
        double r2 = rn.y();
        double r3 = rn.z();
        double r11 = r1 * r1;
        double r22 = r2 * r2;
        double r33 = r3 * r3;
        double norm_rv2 = r11 + r22 + r33;
        double norm_rv = sqrt(norm_rv2);
        if (norm_rv <= 1e-5) {
            v.setZero();
            J.setZero();
            J.rightCols<3>().setIdentity();
            return;
        }
        double r12 = r1 * r2;
        double r13 = r1 * r3;
        double r23 = r2 * r3;
        double atan_tmp = std::atan(norm_rv / rnw);
        double tmp_diag = atan_tmp / norm_rv;
        double tmp = (rnw - tmp_diag) / norm_rv2;
        J << -r1, r11 * tmp + tmp_diag, r12 * tmp, r13 * tmp,
             -r2, r12 * tmp, r22 * tmp + tmp_diag, r23 * tmp,
             -r3, r13 * tmp, r23 * tmp, r33 * tmp + tmp_diag;
        v = tmp_diag * rn.vec();
    }

    static void drot(const Eigen::Vector3d& v,const Eigen::Quaterniond& q, Eigen::Matrix<double, 3, 4> &J)
    {
        double qw = q.w();
        double qx = q.x();
        double qy = q.y();
        double qz = q.z();
        double v1 = v(0);
        double v2 = v(1);
        double v3 = v(2);
        Eigen::Vector3d vec;
        vec << (v1 * qw + v2 * qz - v3 * qy) * 2,
               (v2 * qw - v1 * qz + v3 * qx) * 2,
               (v3 * qw + v1 * qy - v2 * qx) * 2;
        double tmp = (v1 * qx + v2 * qy + v3 * qz) * 2;
        J << vec(0), tmp, -vec(2), vec(1),
             vec(1), vec(2), tmp, -vec(0),
             vec(2), -vec(1), vec(0), tmp;
    }

    static void drot2(const Eigen::Vector3d& v,const Eigen::Quaterniond& q, Eigen::Matrix<double, 3, 4> &J)
    {
        double qw = q.w();
        double qx = q.x();
        double qy = q.y();
        double qz = q.z();
        double v1 = v(0);
        double v2 = v(1);
        double v3 = v(2);
        Eigen::Vector3d vec;
        vec << (v1 * qw - v2 * qz + v3 * qy) * 2,
               (v2 * qw + v1 * qz - v3 * qx) * 2,
               (v3 * qw - v1 * qy + v2 * qx) * 2;
        double tmp = (v1 * qx + v2 * qy + v3 * qz)*2;
        J << vec(0), tmp, vec(2), -vec(1),
             vec(1), -vec(2), tmp, vec(0),
             vec(2), vec(1), -vec(0), tmp;
    }

    static Eigen::Quaterniond deltaQ(const Eigen::Vector3d& theta)
    {
        return Eigen::Quaterniond(1.0, theta(0)/2, theta(1)/2, theta(2)/2);
    }
};

class Sphere
{
public:
    static Eigen::Matrix<double, 3, 2> TangentBasis(const Eigen::Vector3d& g0)
    {
        Eigen::Vector3d b, c;
        Eigen::Vector3d a = g0.normalized();
        Eigen::Vector3d tmp(0, 0, 1);
        if(a != tmp) {
            b = (tmp - a * (a(2))).normalized();
        } else {
            tmp << 1, 0, 0;
            b = (tmp - a * (a(0))).normalized();
        }
        c = a.cross(b);
        Eigen::Matrix<double, 3, 2> bc;
        bc.block<3, 1>(0, 0) = b;
        bc.block<3, 1>(0, 1) = c;
        return bc;
    }
};
