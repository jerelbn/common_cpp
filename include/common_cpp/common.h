// Name: Jerel Nielsen
// Date: 15 June 2017
// Desc: Container for common functions.

#ifndef COMMON_H
#define COMMON_H


#include <iostream>
#include <yaml-cpp/yaml.h>
#include <eigen3/Eigen/Eigen>
#include <random>

namespace common
{


template<typename T>
static const T gravity = 9.80665;

static const Eigen::Vector3d e1 = [] {
  Eigen::Vector3d tmp;
  tmp << 1, 0, 0;
  return tmp;
}();

static const Eigen::Vector3d e2 = [] {
  Eigen::Vector3d tmp;
  tmp << 0, 1, 0;
  return tmp;
}();

static const Eigen::Vector3d e3 = [] {
  Eigen::Vector3d tmp;
  tmp << 0, 0, 1;
  return tmp;
}();

static const Eigen::Matrix3d I_3x3 = [] {
  Eigen::Matrix3d tmp;
  tmp << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  return tmp;
}();

static const Eigen::Matrix2d I_2x2 = [] {
  Eigen::Matrix2d tmp;
  tmp << 1, 0, 0, 1;
  return tmp;
}();

// rotation from NED style camera body coordinates to camera coordinates
static const Eigen::Matrix2d R_cb2c = [] {
  Eigen::Matrix2d tmp;
  tmp << 0, 1, 0,
         0, 0, 1,
         1, 0, 0;
  return tmp;
}();

class Quaternion
{

public:

  Quaternion();
  ~Quaternion();
  Quaternion(double _w, double _x, double _y, double _z);
  Quaternion(double roll, double pitch, double yaw);
  Quaternion(Eigen::Vector4d v);
  Quaternion(Eigen::Vector3d fz);

  Quaternion operator*(const Quaternion &q2);
  Quaternion operator+(const Eigen::Vector3d &delta);
  Eigen::Vector3d operator-(Quaternion &q2);
  friend std::ostream& operator<<(std::ostream &os, const Quaternion &q);

  Quaternion inv();
  double mag();
  void normalize();
  double roll();
  double pitch();
  double yaw();
  void fromEigen(const Eigen::Vector4d q);
  Eigen::Vector4d toEigen();
  Eigen::Vector3d bar();
  Eigen::Matrix3d R();
  Eigen::Vector3d rotSlow(Eigen::Vector3d v);
  Eigen::Vector3d rot(Eigen::Vector3d v);
  Eigen::Vector3d uvec();
  Eigen::MatrixXd proj();
  Quaternion exp(const Eigen::Vector3d delta);
  Eigen::Vector3d log(const Quaternion q);

private:

  double w;
  double x;
  double y;
  double z;

};


// skew symmetric matrix from vector
template<typename T>
Eigen::Matrix<T,3,3> skew(const Eigen::Matrix<T,3,1>& vec)
{
  Eigen::Matrix<T,3,3> A;
  A <<     0, -vec(2),  vec(1),
      vec(2),       0, -vec(0),
     -vec(1),  vec(0),       0;
  return A;
}


// vector from skew symmetric matrix
template<typename T>
Eigen::Matrix<T,3,1> vex(const Eigen::Matrix<T,3,3>& mat)
{
  Eigen::Matrix<T,3,1> v(mat(2,1), mat(0,2), mat(1,0));
  return v;
}


// matrix exponential given skew symmetric delta
template<typename T>
Eigen::Matrix<T,3,3> expR(const Eigen::Matrix<T,3,3>& deltax)
{
  Eigen::Matrix<T,3,1> axis_angle = vex(deltax);
  T theta = axis_angle.norm();
  if (theta > 1e-6)
    return I_3x3 + sin(theta)/theta*deltax + (1-cos(theta))/theta/theta*deltax*deltax;
  else
    return I_3x3;
}


// rotation matrix logarithmic map to vector
template<typename T>
Eigen::Matrix<T,3,3> logR(const Eigen::Matrix<T,3,3>& R)
{
  // rotation magnitude
  T theta = acos((R.trace()-1)/2.0);

  // avoid numerical error with approximation
  Eigen::Matrix<T,3,3> deltax;
  if (theta > 1e-6)
    deltax = theta/(2*sin(theta))*(R - R.transpose());
  else
    deltax = 0.5*(R - R.transpose());

  return deltax;
}


// rotation from vehicle-2 to body frame
template<typename T>
Eigen::Matrix<T,3,3> R_v2_to_b(const T& phi)
{
  Eigen::Matrix<T,3,3> R_v22b;
  R_v22b << 1,         0,        0,
            0,  cos(phi), sin(phi),
            0, -sin(phi), cos(phi);
  return R_v22b;
}

// rotation from vehicle-1 to vehicle-2 frame
template<typename T>
Eigen::Matrix<T,3,3> R_v1_to_v2(const T& theta)
{
  Eigen::Matrix<T,3,3> R_v12v2;
  R_v12v2 << cos(theta), 0, -sin(theta),
                      0, 1,           0,
             sin(theta), 0,  cos(theta);
  return R_v12v2;
}


// rotation from vehicle to vehicle-1 frame
template<typename T>
Eigen::Matrix<T,3,3> R_v_to_v1(const T& psi)
{
  Eigen::Matrix<T,3,3> R_v2v1;
  R_v2v1 <<  cos(psi), sin(psi), 0,
            -sin(psi), cos(psi), 0,
                    0,        0, 1;
  return R_v2v1;
}

// rotation from vehicle to body frame (3-2-1 Euler)
template<typename T>
Eigen::Matrix<T,3,3> R_v_to_b(const T& phi, const T& theta, const T& psi)
{
  return R_v2_to_b(phi) * R_v1_to_v2(theta) * R_v_to_v1(psi);
}


// Loads scalar parameters from a .yaml file
// Author: James Jackson
template <typename T>
bool get_yaml_node(const std::string key, const std::string filename, T& val, bool print_error = true) 
{
  YAML::Node node = YAML::LoadFile(filename);
  if (node[key])
  {
    val = node[key].as<T>();
    return true;
  }
  else
  {
    if (print_error)
    {
      printf("Unable to load \"%s\" from %s\n", key.c_str(), filename.c_str());
    }
    return false;
  }
}


// Loads array from a .yaml file into an Eigen-type matrix or vector.
// Author: James Jackson
template <typename T>
bool get_yaml_eigen(const std::string key, const std::string filename, Eigen::MatrixBase<T>& val) 
{
  YAML::Node node = YAML::LoadFile(filename);
  std::vector<double> vec;
  if (node[key])
  {
    vec = node[key].as<std::vector<double>>();
    if (vec.size() == (val.rows() * val.cols()))
    {
      int k = 0;
      for (int i = 0; i < val.rows(); i++)
      {
        for (int j = 0; j < val.cols(); j++)
        {
          val(i,j) = vec[k++];
        }
      }
      return true;
    }
    else
    {
      printf("Eigen Matrix Size does not match parameter size for \"%s\" in %s", key.c_str(), filename.c_str());
      return false;
    }
  }
  else
  {
    printf("Unable to load \"%s\" from %s\n", key.c_str(), filename.c_str());
    return false;
  }
}


// Saturates a scalar value.
template <typename T>
T saturate(const T& val, const T& max, const T& min)
{
  if (val > max)
    return max;
  if (val < min)
    return min;
  return val;
}


// Random normal matrix generation
template<typename T1, typename T2>
void randomNormalMatrix(Eigen::MatrixBase<T1>& matrix, std::normal_distribution<T2>& dist, std::default_random_engine& rng)
{
  for (int i = 0; i < matrix.rows(); ++i)
    for (int j = 0; j < matrix.cols(); ++j)
      matrix(i,j) = dist(rng);
}


} // namespace common

#endif // COMMON_H
