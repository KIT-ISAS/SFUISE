#include <ros/ros.h>
#include "../include/Linearizer.h"
#include "../include/Accumulator.h"
#include "../include/utils/tic_toc.h"
#include "std_msgs/Int64.h"
#include "sensor_msgs/Imu.h"
#include "cf_msgs/Tdoa.h"
#include "sensor_msgs/PointCloud.h"
#include "isas_msgs/Anchorlist.h"
#include "isas_msgs/RTLSStick.h"
#include "sfuise_msgs/Calib.h"
#include "sfuise_msgs/Spline.h"
#include "sfuise_msgs/Estimate.h"

class SplineFusion
{

  public:

    SplineFusion(ros::NodeHandle& nh)
    {
        if_anchor_ini = false;
        average_runtime = 0;
        window_count = 0;
        solver_flag = INITIAL;
        readParameters(nh);
        sub_imu = nh.subscribe("/EstimationInterface/imu_ds", 1000, &SplineFusion::getImuCallback, this);
        sub_anchor = nh.subscribe("/EstimationInterface/anchor_list", 1000, &SplineFusion::getAnchorCallback, this);
        if (if_tdoa) {
            sub_uwb = nh.subscribe("/EstimationInterface/tdoa_ds", 1000, &SplineFusion::getTdoaCallback, this);
        } else {
            sub_uwb = nh.subscribe("/EstimationInterface/toa_ds", 1000, &SplineFusion::getToaCallback, this);
        }
        pub_knots_active = nh.advertise<sensor_msgs::PointCloud>("active_control_points", 1000);
        pub_knots_inactive = nh.advertise<sensor_msgs::PointCloud>("inactive_control_points", 1000);
        pub_calib = nh.advertise<sfuise_msgs::Calib>("sys_calib", 100);
        pub_est = nh.advertise<sfuise_msgs::Estimate>("est_window", 1000);
        pub_start_time = nh.advertise<std_msgs::Int64>("start_time", 1000);
    }

    void run()
    {
        static int num_window = 0;
        TicToc t_window;
        if (initialization()) {
            displayControlPoints();
            optimization();
            double t_consum = t_window.toc();
            average_runtime = (t_consum + double(num_window) * average_runtime) / double (num_window + 1);
            num_window++;
            if ((int) window_count <= n_window_calib) {
                sfuise_msgs::Calib calib_msg;
                calib_msg.q_nav_uwb.w = calib_param.q_nav_uwb.w();
                calib_msg.q_nav_uwb.x = calib_param.q_nav_uwb.x();
                calib_msg.q_nav_uwb.y = calib_param.q_nav_uwb.y();
                calib_msg.q_nav_uwb.z = calib_param.q_nav_uwb.z();
                calib_msg.t_nav_uwb.x = calib_param.t_nav_uwb[0];
                calib_msg.t_nav_uwb.y = calib_param.t_nav_uwb[1];
                calib_msg.t_nav_uwb.z = calib_param.t_nav_uwb[2];
                geometry_msgs::Point offset_msg;
                offset_msg.x = calib_param.offset.x();
                offset_msg.y = calib_param.offset.y();
                offset_msg.z = calib_param.offset.z();
                calib_msg.t_tag_body_set = offset_msg;
                pub_calib.publish(calib_msg);
            }
            if (spline_local.numKnots() >= (size_t) window_size) {
                window_count++;
                if (solver_flag == INITIAL) {
                    solver_flag = FULLSIZE;
                }
            }
            sfuise_msgs::Spline spline_msg;
            spline_local.getSplineMsg(spline_msg);
            sfuise_msgs::Estimate est_msg;
            est_msg.spline = spline_msg;
            est_msg.if_full_window.data = (solver_flag != INITIAL);
            est_msg.runtime.data = average_runtime;
            pub_est.publish(est_msg);
            displayControlPoints();
            if (solver_flag == FULLSIZE) spline_local.removeOneOldState();
        }
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  private:

    static constexpr double NS_TO_S = 1e-9;

    ros::Subscriber sub_imu;
    ros::Subscriber sub_anchor;
    ros::Subscriber sub_uwb;
    ros::Publisher pub_knots_active;
    ros::Publisher pub_knots_inactive;
    ros::Publisher pub_calib;
    ros::Publisher pub_est;
    ros::Publisher pub_start_time;

    Parameters param;
    CalibParam calib_param;

    Eigen::aligned_deque<TOAData> toa_buff;
    Eigen::aligned_deque<TDOAData> tdoa_buff;
    Eigen::aligned_deque<ImuData> imu_buff;
    Eigen::aligned_deque<ImuData> imu_window;
    Eigen::aligned_deque<TOAData> toa_window;
    Eigen::aligned_deque<TDOAData> tdoa_window;

    bool if_anchor_ini;
    bool if_tdoa;
    bool if_uwb_only;

    size_t window_count;
    int window_size;
    int n_window_calib;

    int64_t dt_ns;
    int64_t bag_start_time;
    int64_t last_imu_t_ns;
    int64_t next_knot_TimeNs;

    enum SolverFlag {
        INITIAL,
        FULLSIZE
    };
    SolverFlag solver_flag;
    SplineState spline_local;

    size_t bias_block_offset;
    size_t gravity_block_offset;
    size_t hess_size;
    bool pose_fixed;
    int max_iter;
    double lambda;
    double lambda_vee;
    double average_runtime;
    std::vector<double> v_toa_offset;

    void readParameters(ros::NodeHandle& nh)
    {
        if (CommonUtils::readParam<double>(nh, "imu_sample_coeff")==0) {
            if_uwb_only = true;
        } else {
            if_uwb_only = false;
        }
        param.if_opt_g = true;
        param.if_opt_transform = true;
        param.w_uwb = CommonUtils::readParam<double>(nh, "w_uwb");
        max_iter = CommonUtils::readParam<int>(nh, "max_iter");
        dt_ns = 1e9 / CommonUtils::readParam<int>(nh, "control_point_fps");
        if_tdoa = CommonUtils::readParam<bool>(nh, "if_tdoa");
        bag_start_time = 0;
        n_window_calib = CommonUtils::readParam<int>(nh, "n_window_calib");
        window_size = CommonUtils::readParam<int>(nh, "window_size");
        if (n_window_calib == 0) {
            ROS_ERROR_STREAM("n_window_calib cannot be set 0.");
            exit(1);
        } else {
            param.q_nav_uwb_init.setIdentity();
            param.t_nav_uwb_init.setZero();
        }
        std::vector<double> accel_var_inv = CommonUtils::readParam<std::vector<double>>(nh, "accel_var_inv");
        param.accel_var_inv << accel_var_inv.at(0), accel_var_inv.at(1), accel_var_inv.at(2);
        std::vector<double> bias_accel_var_inv = CommonUtils::readParam<std::vector<double>>(nh, "bias_accel_var_inv");
        param.bias_accel_var_inv << bias_accel_var_inv.at(0), bias_accel_var_inv.at(1), bias_accel_var_inv.at(2);
        param.w_acc = CommonUtils::readParam<double>(nh, "w_accel");
        param.w_bias_acc = CommonUtils::readParam<double>(nh, "w_bias_accel");
        std::vector<double> gyro_var_inv = CommonUtils::readParam<std::vector<double>>(nh, "gyro_var_inv");
        param.gyro_var_inv << gyro_var_inv.at(0), gyro_var_inv.at(1), gyro_var_inv.at(2);
        std::vector<double> bias_gyro_var_inv = CommonUtils::readParam<std::vector<double>>(nh, "bias_gyro_var_inv");
        param.bias_gyro_var_inv << bias_gyro_var_inv.at(0), bias_gyro_var_inv.at(1), bias_gyro_var_inv.at(2);
        param.w_gyro = CommonUtils::readParam<double>(nh, "w_gyro");
        param.w_bias_gyro = CommonUtils::readParam<double>(nh, "w_bias_gyro");
        param.if_reject_uwb = CommonUtils::readParam<bool>(nh, "if_reject_uwb");
        if (param.if_reject_uwb) {
            param.reject_uwb_thresh = CommonUtils::readParam<double>(nh, "reject_uwb_thresh");
            param.reject_uwb_window_width = CommonUtils::readParam<double>(nh, "reject_uwb_window_width");
        }
        std::vector<double> v_offset;
        nh.getParam("offset", v_offset);
        calib_param.offset = Eigen::Vector3d(v_offset.at(0), v_offset.at(1), v_offset.at(2));
        if (!if_tdoa) {
            v_toa_offset = CommonUtils::readParam<std::vector<double>>(nh, "toa_offset");
        }
    }

    void getImuCallback(const sensor_msgs::ImuConstPtr& imu_msg)
    {
        int64_t t_ns = imu_msg->header.stamp.toNSec();
        Eigen::Vector3d acc(imu_msg->linear_acceleration.x, imu_msg->linear_acceleration.y, imu_msg->linear_acceleration.z);
        Eigen::Vector3d gyro(imu_msg->angular_velocity.x, imu_msg->angular_velocity.y, imu_msg->angular_velocity.z);
        ImuData imu(t_ns, gyro, acc);
        imu_buff.push_back(imu);
    }

    void getTdoaCallback(const cf_msgs::Tdoa::ConstPtr& msg)
    {
        TDOAData uwb(msg->header.stamp.toNSec(), msg->idA, msg->idB, msg->data);
        tdoa_buff.push_back(uwb);
    }

    void getToaCallback(const isas_msgs::RTLSStick::ConstPtr& uwb_msg)
    {
        int64_t t_ns = uwb_msg->header.stamp.toNSec();
        for (const auto& rg : uwb_msg->ranges) {
            if (rg.ra == 0) continue;
            TOAData uwb(t_ns, rg.id, rg.range);
            toa_buff.push_back(uwb);
        }
    }

    void getAnchorCallback(const isas_msgs::Anchorlist::ConstPtr& anchor_msg)
    {
        if (if_anchor_ini) return;
        for (const auto& anchor : anchor_msg->anchor) {
            param.anchor_list[anchor.id] = Eigen::Vector3d(anchor.position.x, anchor.position.y, anchor.position.z);
        }
        if_anchor_ini = true;
        if (!if_tdoa) {
            int i = 0;
            for (auto it = param.anchor_list.begin(); it != param.anchor_list.end(); it++) {
                param.toa_offset[it->first] = v_toa_offset[i];
                i++;
            }
        }
    }

    template <typename type_data>
    void updateMeasurements(Eigen::aligned_deque<type_data>& data_window, Eigen::aligned_deque<type_data>& data_buff)
    {
        int64_t t_window_l = spline_local.minTimeNs();
        if(!data_window.empty()) {
            while (data_window.front().time_ns < t_window_l) {
               data_window.pop_front();
            }
        }
        int64_t t_window_r = spline_local.maxTimeNs();
        for (size_t i = 0; i < data_buff.size(); i++) {
            auto v = data_buff.at(i);
            if (v.time_ns >= t_window_l && v.time_ns <= t_window_r) {
                data_window.push_back(v);
            } else if (v.time_ns > t_window_r) {
                break;
            }
        }
        while (data_buff.front().time_ns <= t_window_r) {
            data_buff.pop_front();
            if(data_buff.empty()) break;
        }
    }

    bool initialization()
    {
        static bool param_set = false;
        static bool initialize_control_point = false;
        if (initialize_control_point) {
            int64_t min_time = 1e18;
            if (!imu_buff.empty()) min_time = imu_buff.back().time_ns;
            if (!toa_buff.empty()) min_time = std::min(toa_buff.back().time_ns, min_time);
            if (!tdoa_buff.empty()) min_time = std::min(tdoa_buff.back().time_ns, min_time);
            if (min_time > spline_local.nextMaxTimeNs()) {
                Eigen::Quaterniond q_ini = spline_local.getKnotOrt(spline_local.numKnots() - 1);
                Eigen::Quaterniond q_ini_backup = q_ini;
                Eigen::Vector3d pos_ini = spline_local.getKnotPos(spline_local.numKnots() - 1);
                Eigen::Matrix<double, 6, 1> bias_ini = spline_local.getKnotBias(spline_local.numKnots() - 1);
                if (!if_uwb_only) {
                    if (spline_local.numKnots() <=2) {
                        last_imu_t_ns = bag_start_time;
                    } else {
                        last_imu_t_ns = imu_window.back().time_ns;
                    }
                    integration(next_knot_TimeNs, q_ini, pos_ini);
                } else {
                    if (if_tdoa) {
                        pos_ini = tdoaMultilateration(next_knot_TimeNs * NS_TO_S);
                        q_ini  = Eigen::Quaterniond::Identity();
                    } else {
                        ROS_ERROR_STREAM("UWB-only tracking only supported for TDOA data!");
                        exit(1);
                    }
                }
                if (q_ini_backup.dot(q_ini) < 0) q_ini = Eigen::Quaterniond(-q_ini.w(), -q_ini.x(), -q_ini.y(), -q_ini.z());
                spline_local.addOneStateKnot(q_ini, pos_ini, bias_ini);
                next_knot_TimeNs += dt_ns;
                return true;
            } else {
                return false;
            }
        } else {
            if (!param_set) {
                param_set = setParameters();
                std_msgs::Int64 start_time;
                start_time.data = bag_start_time;
                pub_start_time.publish(start_time);
            }
            if (param_set && if_anchor_ini) {
                spline_local.init(dt_ns, 0, bag_start_time);
                if (!if_uwb_only) {
                    Eigen::Vector3d gravity_sum(0, 0, 0);
                    size_t n_imu = imu_buff.size();
                    for (size_t i = 0; i < n_imu; i++) {
                        gravity_sum += imu_buff.at(i).accel;
                    }
                    gravity_sum /= n_imu;
                    Eigen::Vector3d gravity_ave = gravity_sum.normalized() * 9.81;
                    calib_param.gravity = gravity_ave;
                }
                calib_param.q_nav_uwb = param.q_nav_uwb_init;
                calib_param.t_nav_uwb = param.t_nav_uwb_init;
                initialize_control_point = true;
                int num = 1;
                for (int i = 0; i < num; i++) {
                    Eigen::Quaterniond q_ini = Eigen::Quaterniond::Identity();
                    Eigen::Vector3d pos_ini = Eigen::Vector3d::Zero();
                    Eigen::Matrix<double, 6, 1> bias_ini = Eigen::Matrix<double, 6, 1>::Zero();
                    spline_local.addOneStateKnot(q_ini, pos_ini, bias_ini);
                    next_knot_TimeNs += dt_ns;
                }
            }
            return false;
        }
    }

    Eigen::Vector3d tdoaMultilateration(double t_s) const
    {
        static Eigen::Vector3d last_knot(0, 0, 0);
        static bool set_origin = true;
        int num_data = 7;
        size_t idx[num_data];
        findClosestNWithOrderedID(t_s, num_data, idx);
        Eigen::Vector3d pos;
        Eigen::MatrixXd H(num_data, 4);
        Eigen::VectorXd b(num_data);
        Eigen::Vector3d anchor0 = param.anchor_list.at(0);
        std::vector<double> range;
        range.push_back(tdoa_buff.at(idx[0]).data);
        for (int i = 1; i < num_data; i++) {
            range.push_back(range[i - 1] + tdoa_buff.at(idx[i]).data);
        }
        for (int i = 0; i < num_data; i++) {
            TDOAData uwb = tdoa_buff.at(idx[i]);
            double range0 = range[i];
            Eigen::Vector3d anchori = param.anchor_list.at(uint16_t(uwb.idB));
            H.row(i) = Eigen::Vector4d(anchori.x() - anchor0.x(), anchori.y() - anchor0.y(), anchori.z() - anchor0.z(), range0);
            b[i] = range0 * range0 - anchori.dot(anchori) + anchor0.dot(anchor0);
        }
        H *= -2;
        Eigen::VectorXd x = (H.transpose() * H).inverse() * H.transpose() * b;
        pos =  x.head(3);
        if (!set_origin) {
            pos.setZero();
            set_origin = true;
        } else {
            pos = calib_param.q_nav_uwb.inverse() * (pos - calib_param.t_nav_uwb);
        }
        return pos;
    }

    void findClosestNWithOrderedID(double t_s, int N, size_t* idx) const
    {
        std::vector<std::pair<std::pair<int,int>,double>> diff;
        for (auto it = tdoa_buff.begin(); it != tdoa_buff.end(); it++) {
            double t = it->time_ns * NS_TO_S;
            diff.push_back(std::make_pair(std::make_pair(it->idA, it->idB), std::abs(t - t_s)));
        }
        int idA = 0;
        int idB = 1;
        for (int i = 0; i < N; i++) {
            size_t idx_min = 0;
            std::pair<std::pair<int,int>,double> it_min = diff[idx_min];
            for (size_t j = 1; j < diff.size(); j++) {
                std::pair<std::pair<int,int>,double> it = diff[j];
                if (it.first == std::make_pair(idA, idB)) {
                    idx_min = it.second < it_min.second ? j : idx_min;
                    it_min = diff[idx_min];
                }
            }
            idx[i] = idx_min;
            idA++;
            idB++;
            diff[idx_min].second = std::numeric_limits<double>::max();
        }
    }

    bool optimization()
    {
        if (spline_local.numKnots() < 2) return false;
        if (param.if_reject_uwb) {
            if (solver_flag == INITIAL &&
                spline_local.numKnots() < int (param.reject_uwb_window_width * window_size))
                param.if_reject_uwb_in_optimization = false;
            else
                param.if_reject_uwb_in_optimization = true;
        } else {
            param.if_reject_uwb_in_optimization = false;
        }
        if ((int) window_count >= n_window_calib) {
            param.if_opt_g = false;
            param.if_opt_transform = false;
        }
        if (!if_uwb_only) updateMeasurements(imu_window, imu_buff);
        if (if_tdoa) {
            updateMeasurements(tdoa_window, tdoa_buff);
        } else {
            updateMeasurements(toa_window, toa_buff);
        }
        bool converged = false;
        int opt_iter = 0;
        pose_fixed = false;
        if (solver_flag == INITIAL) {
            pose_fixed = true;
        }
        lambda = 1e-6;
        lambda_vee = 2;
        updateLinearizerSize();
        while (!converged && opt_iter < max_iter) {
            converged = optimize(opt_iter);
            opt_iter++;
        }
        return converged;
    }

    bool setParameters()
    {
        if (imu_buff.empty() && !if_uwb_only) return false;
        if (imu_buff.empty() && toa_buff.empty() && tdoa_buff.empty()) {
            return false;
        } else {
            if (!imu_buff.empty()) {
                bag_start_time += imu_buff.front().time_ns;
            } else if (!toa_buff.empty()) {
                bag_start_time += toa_buff.front().time_ns;
            } else {
                bag_start_time += tdoa_buff.front().time_ns;
            }
        }
        next_knot_TimeNs = bag_start_time;
        return true;
    }

    void integrateStep(int64_t prevTime, int64_t dt_, const ImuData& imu, Eigen::Matrix3d& Rs, Eigen::Vector3d& Ps, Eigen::Vector3d& Vs)
    {
        static bool first_imu = false;
        static Eigen::Vector3d acc_0;
        static Eigen::Vector3d gyr_0;
        Eigen::Vector3d linear_acceleration = imu.accel;
        Eigen::Vector3d angular_velocity = imu.gyro;
        if (!first_imu) {
            first_imu = true;
            acc_0 = linear_acceleration;
            gyr_0 = angular_velocity;
        }
        Eigen::Vector3d g = calib_param.gravity;
        Eigen::Vector3d ba;
        Eigen::Vector3d bg;
        Eigen::Matrix<double, 6, 1> bias = spline_local.itpBias(prevTime);
        ba = bias.head<3>();
        bg = bias.tail<3>();
        double dt = dt_ * NS_TO_S;
        Eigen::Vector3d un_acc_0;
        un_acc_0 = Rs * (acc_0 - ba) - g;
        Eigen::Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - bg;
        Rs *= Quater::deltaQ(un_gyr * dt).toRotationMatrix();
        Eigen::Vector3d un_acc_1 = Rs * (linear_acceleration - ba) - g;
        Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
        Ps += dt * Vs + 0.5 * dt * dt * un_acc;
        Vs += dt * un_acc;
        acc_0 = linear_acceleration;
        gyr_0 = angular_velocity;
    }

    void integration(const int64_t curTime, Eigen::Quaterniond& qs, Eigen::Vector3d& Ps)
    {
        std::vector<ImuData> imu_vec;
        getIMUInterval(last_imu_t_ns, curTime, imu_vec);
        if (!imu_vec.empty()) {
            Eigen::Quaterniond qs0;
            spline_local.itpQuaternion(last_imu_t_ns, &qs0);
            Eigen::Matrix3d Rs0(qs0);
            Eigen::Vector3d Ps0 = spline_local.itpPosition(last_imu_t_ns);
            Eigen::Vector3d Vs0 = spline_local.itpPosition<1>(last_imu_t_ns);
            for(size_t i = 0; i < imu_vec.size(); i++) {
                int64_t dt;
                int64_t t_ns = imu_vec[i].time_ns;
                if(i == 0) {
                    dt = t_ns - last_imu_t_ns;
                } else {
                    dt = t_ns - imu_vec[i - 1].time_ns;
                }
                integrateStep(last_imu_t_ns, dt, imu_vec[i], Rs0, Ps0, Vs0);
            }
            qs = Eigen::Quaterniond(Rs0);
            Ps = Ps0;
        } else {
            qs = spline_local.extrapolateOrtKnot(1);
            Ps = spline_local.extrapolatePosKnot(1);
        }
    }

    bool getIMUInterval(int64_t t0, int64_t t1, std::vector<ImuData>& imu_vec)
    {
        if(imu_buff.empty()) {
            printf("No IMU available. \n");
            return false;
        }
        int idx = 0;
        while(imu_buff.at(idx).time_ns <= std::min(imu_buff.back().time_ns, t1)) {
            imu_vec.push_back(imu_buff.at(idx));
            idx++;
            if(idx >= imu_buff.size()) break;
        }
        return true;
    }

    void displayControlPoints()
    {
        sensor_msgs::PointCloud points_inactive_msg;
        points_inactive_msg.header.frame_id = "map";
        points_inactive_msg.header.stamp.fromNSec(spline_local.minTimeNs());
        points_inactive_msg.points.push_back(getPointMsg(spline_local.getIdlePos(0)));
        points_inactive_msg.points.push_back(getPointMsg(spline_local.getIdlePos(1)));
        points_inactive_msg.points.push_back(getPointMsg(spline_local.getIdlePos(2)));
        sensor_msgs::PointCloud points_active_msg;
        points_active_msg.header.frame_id = "map";
        points_active_msg.header.stamp.fromNSec(spline_local.minTimeNs());
        for (size_t i = 0; i < spline_local.numKnots(); i++) {
            points_active_msg.points.push_back(getPointMsg(spline_local.getKnotPos(i)));
        }
        pub_knots_inactive.publish(points_inactive_msg);
        pub_knots_active.publish(points_active_msg);
    }

    geometry_msgs::Point32 getPointMsg(Eigen::Vector3d p)
    {
        geometry_msgs::Point32 p_msg;
        Eigen::Vector3d p_U = calib_param.q_nav_uwb * p + calib_param.t_nav_uwb;
        p_msg.x = p_U.x();
        p_msg.y = p_U.y();
        p_msg.z = p_U.z();
        return p_msg;
    }

    bool optimize(const int iter)
    {
        Linearizer lopt(bias_block_offset, gravity_block_offset,
                        hess_size, &spline_local, &calib_param, &param, pose_fixed);
        if (!imu_window.empty()) lopt(imu_window);
        if (!tdoa_window.empty()) lopt(tdoa_window);
        if (!toa_window.empty()) lopt(toa_window);
        if (iter) {
            double gradient_max_norm = lopt.accum.getB().array().abs().maxCoeff();
            if (gradient_max_norm < 1e-8) return true;
        }
        lopt.accum.setup_solver();
        Eigen::VectorXd Hdiag = lopt.accum.Hdiagonal();
        bool stop = false;
        while (!stop) {
            Eigen::VectorXd Hdiag_lambda = Hdiag * lambda;
            for (int i = 0; i < Hdiag_lambda.size(); i++) {
              Hdiag_lambda[i] = std::max(Hdiag_lambda[i], 1e-18);
            }
            Eigen::VectorXd inc_full = -lopt.accum.solve(&Hdiag_lambda);
            Eigen::aligned_deque<Eigen::Vector3d> knots_trans_backup;
            Eigen::aligned_deque<Eigen::Quaterniond> knots_rot_backup;
            Eigen::aligned_deque<Eigen::Matrix<double, 6, 1>> knots_bias_backup;
            spline_local.getAllStateKnots(knots_trans_backup, knots_rot_backup, knots_bias_backup);
            CalibParam calib_param_backup = calib_param;
            applyIncFull(inc_full);
            ComputeErrorSplineOpt eopt(&spline_local, &calib_param, &param);
            if (!toa_window.empty()) eopt(toa_window);
            if (!imu_window.empty()) eopt(imu_window);
            if (!tdoa_window.empty()) eopt(tdoa_window);
            double f_diff = lopt.error - eopt.error;
            double l_diff = 0.5 * inc_full.dot(inc_full *lambda - lopt.accum.getB());
            double step_quality = f_diff / l_diff;
            if (step_quality < 0) {
                lambda = std::min(100.0, lambda_vee * lambda);
                if (abs(lambda - 100.0) < 1e-3) {
                    stop = true;
                }
                lambda_vee *= 2;
                spline_local.setAllKnots(knots_trans_backup, knots_rot_backup, knots_bias_backup);
                calib_param.setCalibParam(calib_param_backup);
            } else {
                if (inc_full.norm()/((double)spline_local.numKnots()) < 1e-10 || abs(f_diff)/lopt.error < 1e-6) {
                    stop = true;
                }
                lambda = std::max(1e-18, lambda * std::max(1.0 / 3, 1 - std::pow(2 * step_quality - 1, 3.0)));
                lambda_vee = 2;
                break;
            }
        }
        return stop;
    }

    void updateLinearizerSize()
    {
        int num_knots = spline_local.numKnots();
        bias_block_offset = Linearizer::POSE_SIZE * num_knots;
        hess_size = bias_block_offset;
        if (!if_uwb_only) {
            hess_size += Linearizer::ACCEL_BIAS_SIZE * num_knots;
            hess_size += Linearizer::GYRO_BIAS_SIZE * num_knots;
        }
        gravity_block_offset = hess_size;
        hess_size += Linearizer::G_SIZE;
        if (param.if_opt_transform) {
          hess_size += Linearizer::OFFSET_SIZE;
          hess_size += Linearizer::ROTATION_SIZE;
        }
    }

    void applyIncFull(Eigen::VectorXd& inc_full)
    {
        size_t num_knots = spline_local.numKnots();
        for (size_t i = 0; i < num_knots; i++) {
          Eigen::Matrix<double, 6, 1> inc = inc_full.segment<Linearizer::POSE_SIZE>(Linearizer::POSE_SIZE * i);
          spline_local.applyPoseInc(i, inc);
        }
        spline_local.checkQuaternionControlPoints();
        if (!if_uwb_only) {
            for (size_t i = 0; i < num_knots; i++) {
                Eigen::Matrix<double, 6, 1> inc = inc_full.segment<Linearizer::BIAS_SIZE>(bias_block_offset + Linearizer::BIAS_SIZE * i);
                spline_local.applyBiasInc(i, inc);
            }
            spline_local.updateBiasIdleFirstWindow();
            if (param.if_opt_g) {
                Eigen::VectorXd dg = inc_full.segment<Linearizer::G_SIZE>(gravity_block_offset);
                Eigen::Vector3d g0 = (calib_param.gravity + Sphere::TangentBasis(calib_param.gravity) * dg).normalized() * 9.81;
                calib_param.gravity = g0;
            }
        }
        if (param.if_opt_transform) {
            calib_param.t_nav_uwb += inc_full.segment<Linearizer::OFFSET_SIZE>(gravity_block_offset + Linearizer::TRANS_OFFSET);
            Eigen::Quaterniond q_inc;
            Quater::exp(inc_full.segment<Linearizer::ROTATION_SIZE>(gravity_block_offset + Linearizer::ROTATION_OFFSET), q_inc);
            calib_param.q_nav_uwb *= q_inc;
        }
    }

};

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "sfuise");
    ROS_INFO("\033[1;32m---->\033[0m Starting SplineFusion.");
    ros::NodeHandle nh("~");
    SplineFusion estimator(nh);
    ros::Rate rate(1000);
    while (ros::ok()) {
        ros::spinOnce();
        estimator.run();
        rate.sleep();
    }
    return 0;
}
