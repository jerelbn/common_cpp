// Name: Jerel Nielsen
// Date: 15 June 2017
// Desc: Container for common functions.

#ifndef COMMON_H
#define COMMON_H


#include <iostream>
#include <yaml-cpp/yaml.h>
#include <eigen3/Eigen/Eigen>
#include <random>
#include <chrono>

namespace common
{


static const double gravity = 9.80665;

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

static const Eigen::Matrix<double, 2, 3> I_2x3 = [] {
  Eigen::Matrix<double, 2, 3> tmp;
  tmp << 1, 0, 0, 0, 1, 0;
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
  Quaternion(const Eigen::Vector4d &v);
  Quaternion(Eigen::Vector3d &fz);
  Quaternion(std::normal_distribution<double> &dist, std::default_random_engine &rng);

  Quaternion operator*(const Quaternion &q2) const;
  Quaternion operator+(const Eigen::Vector3d &delta) const;
  void operator+=(const Eigen::Vector3d &delta);
  Eigen::Vector3d operator-(const Quaternion &q2) const;
  friend std::ostream& operator<<(std::ostream &os, const Quaternion &q);

  void normalize();
  void fromEigen(const Eigen::Vector4d &q);
  Quaternion inv() const;
  double mag() const;
  double roll() const;
  double pitch() const;
  double yaw() const;
  Eigen::Vector3d euler() const;
  Eigen::Vector4d toEigen() const;
  Eigen::Vector3d bar() const;
  Eigen::Matrix3d R() const;
  Eigen::Vector3d rotSlow(const Eigen::Vector3d &v) const;
  Eigen::Vector3d rot(const Eigen::Vector3d &v) const;
  Eigen::Vector3d uvec() const;
  Eigen::MatrixXd proj() const;
  static Quaternion exp(const Eigen::Vector3d &delta);
  static Eigen::Vector3d log(const Quaternion &q);
  static Eigen::Vector2d log_uvec(const Quaternion &q1, const Quaternion &q2); // q1 - q2
  static Eigen::Matrix3d dexp(const Eigen::Vector3d &delta);

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


// Random normal array generation
template<typename T1, typename T2>
void randomNormalArray(Eigen::ArrayBase<T1>& array, std::normal_distribution<T2>& dist, std::default_random_engine& rng)
{
  for (int i = 0; i < array.rows(); ++i)
    for (int j = 0; j < array.cols(); ++j)
      array(i,j) = dist(rng);
}


// Shows progress of simulation or event
// Author: James Jackson
class ProgressBar
{
public:
  ProgressBar(){}
  ProgressBar(int total, int barwidth) :
    initialized_(false),
    barwidth_(barwidth),
    total_(total)
  {}

  ~ProgressBar()
  {
    std::cout << std::endl;
  }

  void init(int total, int barwidth)
  {
    initialized_ = false;
    barwidth_ = barwidth;
    total_ = total;
    last_completed_ = 0;
  }

  void print(int completed)
  {
    if (!initialized_)
    {
      last_print_time_ = std::chrono::system_clock::now();
      start_time_ = std::chrono::system_clock::now();
      initialized_ = true;
    }
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time_).count()/1000.0;

    // limit printing to about 60 Hz
    if (std::chrono::duration_cast<std::chrono::milliseconds>(now - last_print_time_).count() > 33
        || completed == total_)
    {
      last_print_time_ = now;
      std::cout << "[";
      int pos = barwidth_ * (completed / (double)total_);
      for (int i = 0; i < barwidth_; ++i)
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
      std::cout << "]  ";
      printf("%.0f%% ", (completed / (double)total_)*100.0);
      double it_s = completed / elapsed;
      double left = (total_ - completed) / it_s;
      printf("[%.2f<%.2f, %.2fit/s] \r", elapsed, left, it_s);
    }
    last_completed_ = completed;
  }
private:
  int barwidth_;
  int total_;
  bool initialized_;
  int last_completed_;

  std::chrono::system_clock::time_point start_time_;
  std::chrono::system_clock::time_point last_print_time_;
};


// Perspective projection into an image
template<typename T>
void projToImg(Eigen::Matrix<T,2,1>& pix, const Eigen::Matrix<T,3,1> &lc, const Eigen::Matrix<T,3,3> &K)
{
  pix = K.topRows(2) * (lc / lc(2));
};


// Direction vector from pixel coordinate
template<typename T>
void dirFromPix(Eigen::Matrix<T,3,1> &dir, const Eigen::Matrix<T,2,1> &pix, const Eigen::Matrix<T,3,3> &K_inv)
{
  dir = K_inv * Eigen::Matrix<T,3,1>(pix(0), pix(1), 1);
  dir /= dir.norm();
}


// Angular difference between two vectors
template<typename T>
double angDiffBetweenVecs(const Eigen::Matrix<T,3,1>& v1, const Eigen::Matrix<T,3,1>& v2)
{
  double val = (v1.transpose() * v2)(0) / (v1.norm() * v2.norm());
  if (val > 1)
    return 0;
  else if (val < -1)
    return M_PI;
  else
    return acos(val);
}


// create string of red or green text
std::string redText(const std::string &input);
std::string greenText(const std::string &input);


// test equivalence of Eigen matrices and vectors
template<typename T1, typename T2>
bool TEST(const std::string &test_name, const double &tol, const Eigen::MatrixBase<T1> &mat1, const Eigen::MatrixBase<T2> &mat2)
{
  if (mat1.norm() < 1e-6 || mat2.norm() < 1e-6)
    std::cout << redText("WARNING: ") << "Test values near zero in " << test_name << ".\n";
  if ((mat1 - mat2).norm() < tol)
  {
    std::cout << greenText("[PASSED] ") << test_name << std::endl;
    return true;
  }
  else
  {
    std::cout << redText("[FAILED] ") << test_name << "\n\n";
    std::cout << "Error: " << (mat1 - mat2).norm() << std::endl;
    std::cout << "\nFirst input  = \n" << mat1 << std::endl;
    std::cout << "\nSecond input = \n" << mat2 << std::endl;
    return false;
  }
}


} // namespace common

#endif // COMMON_H
