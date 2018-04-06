// Name: Jerel Nielsen
// Date: 15 June 2017
// Desc: Container for common functions.

#ifndef COMMON_H
#define COMMON_H


#include <iostream>
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

} // namespace common

#endif // COMMON_H
