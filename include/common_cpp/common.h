// Name: Jerel Nielsen
// Date: 15 June 2017
// Desc: Container for common functions.

#ifndef COMMON_H
#define COMMON_H


#include <iostream>
#include <ros/ros.h>
#include <tf/transform_listener.h>
#include <eigen3/Eigen/Eigen>
#include <opencv2/opencv.hpp>

namespace common
{


class Quaternion
{

public:

  Quaternion();
  ~Quaternion();
  Quaternion(double _w, double _x, double _y, double _z);
  Quaternion(double roll, double pitch, double yaw);
  Quaternion(Eigen::Vector4d v);
  Quaternion(Eigen::Vector3d fz);

  double x;
  double y;
  double z;
  double w;

  Quaternion operator*(const Quaternion &q2);
  friend std::ostream& operator<<(std::ostream &os, const Quaternion &q);
  Quaternion inv();
  double mag();
  void normalize();
  double roll();
  double pitch();
  double yaw();
  void convertFromEigen(Eigen::Vector4d q);
  Eigen::Vector4d convertToEigen();
  Eigen::Matrix3d rot();
  Eigen::Vector3d rotateVector(Eigen::Vector3d v);
  Eigen::Vector3d unitVector();
  Eigen::MatrixXd projection();

};

Eigen::VectorXd rk5(Eigen::VectorXd state, Eigen::VectorXd input, std::function<Eigen::VectorXd(Eigen::VectorXd, Eigen::VectorXd)> ode, double h);
Quaternion exp_q(const Eigen::Vector3d delta);
Eigen::Vector3d log_q(const Quaternion q);
Eigen::Vector3d log_R(const Eigen::Matrix3d R);
common::Quaternion vec2quat(const Eigen::Vector3d v);
Eigen::Vector3d vex(const Eigen::Matrix3d mat);
Eigen::Matrix3d skew(const Eigen::Vector3d vec);
Eigen::Matrix3d R_v2_to_b(double phi);
Eigen::Matrix3d R_v1_to_v2(double theta);
Eigen::Matrix3d R_v_to_v1(double psi);
Eigen::Matrix3d R_v_to_b(double phi, double theta, double psi);
Eigen::Matrix3d R_cb2c();
Eigen::MatrixXd cv2eigen(cv::Mat cvmat);
void getQuaternionFromTF(std::string frame1, std::string frame2, common::Quaternion& q);

/*---------------------------------------------
Templated functions implemented in header to avoid 
explicity instantiating all needed template instances.
----------------------------------------------*/

template <typename T1, typename T2>
void rosImportScalar(ros::NodeHandle nh, std::string param, T2& value, T1 default_value)
{
  // get scalar from ROS parameter server
  if (!nh.getParam(param, value))
  {
    ROS_WARN("Could not find %s/%s on the server.", nh.getNamespace().c_str(), param.c_str());
    value = default_value;
  }
}

template <typename T1, typename T2>
void rosImportMatrix(ros::NodeHandle nh, std::string param, Eigen::MatrixBase<T2>& mat)
{
  // get array from ROS parameter server
  std::vector<T1> vec;
  if (!nh.getParam(param, vec))
  {
    ROS_WARN("Could not find %s/%s on the server. Set to zeros.", nh.getNamespace().c_str(), param.c_str());
    mat.setZero();
    return;
  }
 
  // ensure imported array has correct number of values then populate the matrix
  ROS_ASSERT_MSG(vec.size() == mat.rows()*mat.cols(), "Param %s/%s is the wrong size", nh.getNamespace().c_str(),param.c_str());
  for (unsigned i = 0; i < mat.rows(); i++)
    for (unsigned j = 0; j < mat.cols(); j++)
      mat(i,j) = vec[mat.cols()*i+j];
}

template <typename T>
void getRotationFromTF(std::string frame1, std::string frame2, Eigen::MatrixBase<T>& R)
{
  // intantiate a ROS tf listener
  tf::TransformListener tf_listener;

  // lookup the transform (from tf tree) from frame1 to frame2
  tf::StampedTransform transform;
  try 
  {
    tf_listener.waitForTransform(frame1, frame2, ros::Time(0), ros::Duration(10.0) );
    tf_listener.lookupTransform(frame1, frame2, ros::Time(0), transform);
  } 
  catch (tf::TransformException& ex)
  {
    ROS_ERROR("%s",ex.what());
  }

  // convert transform to an Eigen rotation
  tf::Matrix3x3 _R = transform.getBasis();
  for (unsigned i = 0; i < 3; i++)
  {
    for (unsigned j = 0; j < 3; j++)
      R(i,j) = _R[i][j];
  }
}

template <typename T>
void getTranslationFromTF(std::string frame1, std::string frame2, Eigen::MatrixBase<T>& t)
{
  // intantiate a ROS tf listener
  tf::TransformListener tf_listener;

  // lookup the transform (from tf tree) from frame1 to frame2
  tf::StampedTransform transform;
  try 
  {
    tf_listener.waitForTransform(frame1, frame2, ros::Time(0), ros::Duration(10.0) );
    tf_listener.lookupTransform(frame1, frame2, ros::Time(0), transform);
  } 
  catch (tf::TransformException& ex)
  {
    ROS_ERROR("%s",ex.what());
  }

  // convert transform to an Eigen translation from frame1 to frame2 given in frame2
  tf::Vector3 _t = transform.getOrigin();
  for (unsigned i = 0; i < 3; i++)
    t(i) = _t[i];
}

} // namespace common

#endif // COMMON_H
