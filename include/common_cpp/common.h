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

Eigen::VectorXd rk5(Eigen::VectorXd state, Eigen::VectorXd input, std::function<Eigen::VectorXd(Eigen::VectorXd, Eigen::VectorXd)> ode, double h);
Eigen::Matrix3d expR(const Eigen::Matrix3d deltax);
Eigen::Vector3d logR(const Eigen::Matrix3d R);
Eigen::Vector3d vex(const Eigen::Matrix3d mat);
Eigen::Matrix3d skew(const Eigen::Vector3d vec);
Eigen::Matrix3d R_v2_to_b(double phi);
Eigen::Matrix3d R_v1_to_v2(double theta);
Eigen::Matrix3d R_v_to_v1(double psi);
Eigen::Matrix3d R_v_to_b(double phi, double theta, double psi);

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
T saturate(const T val, const T max, const T min)
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
