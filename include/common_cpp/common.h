// Name: Jerel Nielsen
// Date: 15 June 2017
// Desc: Container for common functions.

#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>
#include <eigen3/Eigen/Eigen>
#include <random>
#include <chrono>
#include <cmath>

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
static const Eigen::Matrix3d R_cb2c = [] {
  Eigen::Matrix3d tmp;
  tmp << 0, 1, 0,
         0, 0, 1,
         1, 0, 0;
  return tmp;
}();


// wrap angle to +- input bound (typically [0,2*pi] or [-pi,pi])
template<typename T>
T wrapAngle(const T &angle, const T &bound)
{
  if (angle > bound)
    return angle - T(2.0) * T(M_PI);
  if (angle < bound - T(2.0) * T(M_PI))
    return angle + T(2.0) * T(M_PI);
  return angle;
}


// round to desired decimal place
template<typename T>
T decRound(const T &number, const int &decimal_place)
{
  return round(number * decimal_place) / decimal_place;
}


// skew symmetric matrix from vector
template<typename T>
Eigen::Matrix<T,3,3> skew(const Eigen::Matrix<T,3,1>& vec)
{
  Eigen::Matrix<T,3,3> A;
  A <<    T(0), -vec(2),  vec(1),
        vec(2),    T(0), -vec(0),
       -vec(1),  vec(0),    T(0);
  return A;
}


// vector from skew symmetric matrix
template<typename T>
Eigen::Matrix<T,3,1> vex(const Eigen::Matrix<T,3,3>& mat)
{
  const Eigen::Matrix<T,3,1> v(mat(2,1), mat(0,2), mat(1,0));
  return v;
}


// matrix exponential given skew symmetric delta
template<typename T>
Eigen::Matrix<T,3,3> expR(const Eigen::Matrix<T,3,3>& deltax)
{
  const Eigen::Matrix<T,3,1> axis_angle = vex(deltax);
  const T theta = axis_angle.norm();
  if (theta > 1e-5)
    return I_3x3.cast<T>() + sin(theta) / theta * deltax +
           (T(1.0) - cos(theta)) / theta / theta * deltax * deltax;
  else
    return I_3x3.cast<T>();
}


// rotation matrix logarithmic map
template<typename T>
Eigen::Matrix<T,3,3> logR(const Eigen::Matrix<T,3,3>& R)
{
  // rotation magnitude
  const T theta = acos((R.trace()-1)/2.0);

  // avoid numerical error with approximation
  Eigen::Matrix<T,3,3> deltax;
  if (theta > 1e-5)
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


// copy data to Eigen matrix
template<typename T>
void copy_ptr_to_eigen(const T* ptr, Eigen::Matrix<T,-1,-1>& m)
{
  int i(0), j(0);
  for (int k = 0; k < m.rows() * m.cols(); ++k)
  {
    m(i,j) = ptr[k];
    ++i;
    if (i == m.rows()) { i = 0; ++j; }
  }
}


// Loads array from a binary file onto the heap (in case the file is large)
// and returns the pointer to the beginning of the array
template<typename T>
T* load_binary(const std::string& filename, long& array_size)
{
  std::ifstream file(filename, std::ios::in|std::ios::binary|std::ios::ate);
  if (file.is_open())
  {
    // File opened with ios::ate flag to position the get pointer at end of file
    // then tellg() gets the file size
    std::ifstream::pos_type size = file.tellg();
    char* memblock = new char[size]; // Allocate memory on the heap for the file
    file.seekg(0, std::ios::beg); // Set the get pointer to beginning of file
    file.read(memblock, size); // Read entire file
    file.close();
    std::cout << "Loaded file: " << filename << std::endl;

    // Return array as type T and its size
    array_size = (long)size / sizeof(T);
    return (T*)memblock;
  }
  else
  {
    std::cout << "Unable to open file: " << filename << std::endl;
    exit(0);
  }
}

// Load binary file directly into an Eigen array
template<typename T>
Eigen::Matrix<T,-1,-1> load_binary_to_matrix(const std::string& filename, const int& matrix_rows)
{
  long array_size;
  T* ptr = load_binary<T>(filename, array_size);
  Eigen::Matrix<T,-1,-1> m(matrix_rows, array_size/matrix_rows);
  copy_ptr_to_eigen(ptr, m);
  delete[] ptr;
  return m;
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


// create string of red text
inline std::string redText(const std::string& input)
{
  return "\033[31m" + input + "\033[0m";
}


// create string of green text
inline std::string greenText(const std::string& input)
{
  return "\033[32m" + input + "\033[0m";
}


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
T angDiffBetweenVecs(const Eigen::Matrix<T,3,1>& v1, const Eigen::Matrix<T,3,1>& v2)
{
  T val = (v1.transpose() * v2)(0) / (v1.norm() * v2.norm());
  if (val > 1)
    return 0;
  else if (val < -1)
    return M_PI;
  else
    return acos(val);
}


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
    std::cout << redText("[FAILED] ") << test_name << std::endl;
    std::cout << "\nError: " << (mat1 - mat2).norm() << std::endl;
    std::cout << "\nValue1 = \n" << mat1 << std::endl;
    std::cout << "\nValue2 = \n" << mat2 << std::endl;
    std::cout << "\nDifference = \n" << mat1 - mat2 << std::endl;
    return false;
  }
}


} // namespace common

#endif // COMMON_H
