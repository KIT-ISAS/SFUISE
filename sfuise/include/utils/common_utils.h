#pragma once

#include <ros/ros.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <deque>
#include <map>
#include <unordered_map>
#include <vector>
#include "geometry_msgs/PoseStamped.h"

namespace Eigen {

template <typename T>
using aligned_vector = std::vector<T, Eigen::aligned_allocator<T>>;

template <typename T>
using aligned_deque = std::deque<T, Eigen::aligned_allocator<T>>;

template <typename K, typename V>
using aligned_map = std::map<K, V, std::less<K>,
                    Eigen::aligned_allocator<std::pair<K const, V>>>;

template <typename K, typename V>
using aligned_unordered_map = std::unordered_map<K, V, std::hash<K>, std::equal_to<K>,
                              Eigen::aligned_allocator<std::pair<K const, V>>>;

}

struct PoseData {
  int64_t time_ns;
  Eigen::Quaterniond orient;
  Eigen::Vector3d pos;
  PoseData(){}
  PoseData(int64_t s, Eigen::Quaterniond& q, Eigen::Vector3d& t) : time_ns(s), orient(q), pos(t) {}
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct ImuData {
    const int64_t time_ns;
    const Eigen::Vector3d gyro;
    const Eigen::Vector3d accel;
    ImuData(const int64_t s, const Eigen::Vector3d& w, const Eigen::Vector3d& a)
      : time_ns(s), gyro(w), accel(a) {}
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct TOAData {
    const int64_t time_ns;
    const int anchor_id;
    const double data;
    TOAData(const int64_t s, const int i, const double r) : time_ns(s), anchor_id(i), data(r) {};
};

struct TDOAData {
  const int64_t time_ns;
  const int idA;
  const int idB;
  const double data;
  TDOAData(const int64_t s, const int idxA, const int idxB, const double r) : time_ns(s), idA(idxA), idB(idxB), data(r) {};
};

class CalibParam
{
  public:
    Eigen::Vector3d offset;
    Eigen::Quaterniond q_nav_uwb;
    Eigen::Vector3d t_nav_uwb;
    Eigen::Vector3d gravity;

    CalibParam(){};

    void setCalibParam(CalibParam calib_param)
    {
        offset = calib_param.offset;
        q_nav_uwb = calib_param.q_nav_uwb;
        t_nav_uwb = calib_param.t_nav_uwb;
        gravity =calib_param.gravity;
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class CommonUtils
{
public:
    template <typename T>
    static T readParam(ros::NodeHandle &nh, std::string name)
    {
        T ans;
        if (!nh.getParam(name, ans))
        {
            ROS_ERROR_STREAM("Failed to load " << name);
            nh.shutdown();
        }
        return ans;
    }

    static geometry_msgs::PoseStamped pose2msg(const int64_t t, const Eigen::Vector3d& pos,
                                              const Eigen::Quaterniond& orient)
    {
        geometry_msgs::PoseStamped msg;
        msg.header.stamp.fromNSec(t);
        msg.pose.position.x = pos.x();
        msg.pose.position.y = pos.y();
        msg.pose.position.z = pos.z();
        msg.pose.orientation.w = orient.w();
        msg.pose.orientation.x = orient.x();
        msg.pose.orientation.y = orient.y();
        msg.pose.orientation.z = orient.z();
        return msg;
    }
};
