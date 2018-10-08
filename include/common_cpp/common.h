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
    std::cout << redText("[FAILED] ") << test_name << std::endl;
    std::cout << "\nvalue1 = \n" << mat1 << std::endl;
    std::cout << "\nvalue2 = \n" << mat2 << std::endl;
    return false;
  }
}


template<typename T = double>
class Quaternion
{

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Quaternion()
  {
    arr(0) = T(1.0);
    arr(1) = T(0.0);
    arr(2) = T(0.0);
    arr(3) = T(0.0);
  }

  Quaternion(const T& _w, const T& _x, const T& _y, const T& _z)
  {
    arr(0) = _w;
    arr(1) = _x;
    arr(2) = _y;
    arr(3) = _z;
  }

  Quaternion(const T& roll, const T& pitch, const T& yaw)
  {
    // Pre-calculations
    const T r_2 = roll / T(2.0);
    const T p_2 = pitch / T(2.0);
    const T y_2 = yaw / T(2.0);
    const T sr = sin(r_2);
    const T sp = sin(p_2);
    const T sy = sin(y_2);
    const T cr = cos(r_2);
    const T cp = cos(p_2);
    const T cy = cos(y_2);

    // Output
    arr(0) = cr*cp*cy + sr*sp*sy;
    arr(1) = sr*cp*cy - cr*sp*sy;
    arr(2) = cr*sp*cy + sr*cp*sy;
    arr(3) = cr*cp*sy - sr*sp*cy;
  }

  Quaternion(const Eigen::Matrix<T,4,1>& v)
  {
    arr(0) = v(0);
    arr(1) = v(1);
    arr(2) = v(2);
    arr(3) = v(3);
  }

  // create quaternion from a unit vector
  Quaternion(Eigen::Matrix<T,3,1>& fz)
  {
    // convert to axis-angle representation
    fz.normalize(); // enforce unit length
    const T theta = acos(fz.dot(e3.cast<T>()));
    if (theta < T(1e-6))
    {
      arr(0) = T(1.0);
      arr(1) = T(0.0);
      arr(2) = T(0.0);
      arr(3) = T(0.0);
    }
    else
    {
      const Eigen::Matrix<T,3,1> iaa = (fz.cross(e3.cast<T>())).normalized();
      const T theta_2 = theta / T(2.0);
      const Eigen::Matrix<T,3,1> qv = iaa * sin(theta_2);

      arr(0) = cos(theta_2);
      arr(1) = qv(0);
      arr(2) = qv(1);
      arr(3) = qv(2);
    }
  }

  void normalize()
  {
    const T m = mag();
    arr(0) /= m;
    arr(1) /= m;
    arr(2) /= m;
    arr(3) /= m;
  }

  // initialize random unit quaternion
  Quaternion(std::normal_distribution<T>& dist, std::default_random_engine& rng)
  {
    Quaternion<T> q(dist(rng), dist(rng), dist(rng), dist(rng));
    q.normalize();
    if (q.w() < 0) q.setW(-q.w());
    arr(0) = q.w();
    arr(1) = q.x();
    arr(2) = q.y();
    arr(3) = q.z();
  }

  // overload addition operator as boxplus for a quaternion and a 3-vector
  Quaternion<T> operator+(const Eigen::Matrix<T,3,1>& delta) const
  {
    return *this * exp(delta);
  }

  void operator+=(const Eigen::Matrix<T,3,1>& delta)
  {
    *this = *this + delta;
  }

  // overload minus operator as boxminus for two quaternions
  Eigen::Matrix<T,3,1> operator-(const Quaternion<T> &q2) const
  {
    return log(q2.inv() * *this);
  }

  Quaternion<T> operator*(const Quaternion<T> &q2) const
  {
    const T qw = w()*q2.w() - x()*q2.x() - y()*q2.y() - z()*q2.z();
    const T qx = w()*q2.x() + x()*q2.w() + y()*q2.z() - z()*q2.y();
    const T qy = w()*q2.y() - x()*q2.z() + y()*q2.w() + z()*q2.x();
    const T qz = w()*q2.z() + x()*q2.y() - y()*q2.x() + z()*q2.w();
    return Quaternion<T>(qw, qx, qy, qz);
  }

  void operator*=(const Quaternion<T> &q)
  {
    *this = *this * q;
  }

  friend std::ostream& operator<<(std::ostream &os, const Quaternion<T> &q)
  {
    os << q.w() << "\n" << q.x() << "\n" << q.y() << "\n" << q.z();
    return os;
  }

  template<typename T2>
  Quaternion<T2> cast() const
  {
    return Quaternion<T2>(arr.cast<T2>());
  }

  void scale(const T& s)
  {
    arr(0) *= s;
    arr(1) *= s;
    arr(2) *= s;
    arr(3) *= s;
  }

  T mag() const
  {
    return sqrt(w()*w() + x()*x() + y()*y() + z()*z());
  }

  T roll() const
  {
    return atan2(T(2.0) * (w()*x() + y()*z()), T(1.0) - T(2.0) * (x()*x() + y()*y()));
  }

  T pitch() const
  {
    const T val = T(2.0) * (w()*y() - x()*z());

    // hold at 90 degrees if invalid
    if (fabs(val) > T(1.0))
      return copysign(T(1.0), val) * T(M_PI) / T(2.0);
    else
      return asin(val);
  }

  T yaw() const
  {
    return atan2(T(2.0) * (w()*z() + x()*y()), T(1.0) - T(2.0) * (y()*y() + z()*z()));
  }

  Eigen::Matrix<T,3,1> euler() const
  {
    return Eigen::Matrix<T,3,1>(roll(),pitch(),yaw());
  }

  void fromEigen(const Eigen::Matrix<T,4,1>& q)
  {
    arr(0) = q(0);
    arr(1) = q(1);
    arr(2) = q(2);
    arr(3) = q(3);
  }

  Eigen::Matrix<T,4,1> toEigen() const
  {
    return Eigen::Matrix<T,4,1>(w(),x(),y(),z());
  }

  Eigen::Matrix<T,3,1> bar() const
  {
    return Eigen::Matrix<T,3,1>(x(),y(),z());
  }

  Quaternion<T> inv() const
  {
    return Quaternion<T>(w(), -x(), -y(), -z());
  }

  Eigen::Matrix<T,3,3> R() const
  {
    // Pre-calculations
    const T ww = w() * w();
    const T wx = w() * x();
    const T wy = w() * y();
    const T wz = w() * z();
    const T xy = x() * y();
    const T xz = x() * z();
    const T yz = y() * z();

    // Output
    Eigen::Matrix<T,3,3> R;
    R(0,0) = T(2.0) * (ww + x() * x()) - T(1.0);
    R(0,1) = T(2.0) * (wz + xy);
    R(0,2) = T(2.0) * (xz - wy);
    R(1,0) = T(2.0) * (xy - wz);
    R(1,1) = T(2.0) * (ww + y() * y()) - T(1.0);
    R(1,2) = T(2.0) * (wx + yz);
    R(2,0) = T(2.0) * (wy + xz);
    R(2,1) = T(2.0) * (yz - wx);
    R(2,2) = T(2.0) * (ww + z() * z()) - T(1.0);
    return R;
  }

  Eigen::Matrix<T,3,1> rotSlow(const Eigen::Matrix<T,3,1>& v) const
  {
    const Quaternion qv(T(0.0), v(0), v(1), v(2));
    const Quaternion qv_new = inv() * qv * *this;
    return Eigen::Matrix<T,3,1>(qv_new.x(), qv_new.y(), qv_new.z());
  }

  Eigen::Matrix<T,3,1> rot(const Eigen::Matrix<T,3,1>& v) const
  {
    const Eigen::Matrix<T,3,1> t = T(2.0) * v.cross(bar());
    return v + w() * t + t.cross(bar());
  }

  Eigen::Matrix<T,3,1> uvec() const
  {
    return rot(e3.cast<T>());
  }

  Eigen::Matrix<T,3,2> proj() const
  {
    return R() * I_2x3.cast<T>().transpose();
  }

  static Quaternion<T> exp(const Eigen::Matrix<T,3,1>& delta)
  {
    const T delta_norm = delta.norm();

    Quaternion<T> q;
    if (delta_norm < T(1e-6)) // avoid numerical error with approximation
    {
      q.setW(T(1.0));
      q.setX(delta(0) / T(2.0));
      q.setY(delta(1) / T(2.0));
      q.setZ(delta(2) / T(2.0));
    }
    else
    {
      const T delta_norm_2 = delta_norm / T(2.0);
      const T sn = sin(delta_norm_2) / delta_norm;
      q.setW(cos(delta_norm_2));
      q.setX(sn * delta(0));
      q.setY(sn * delta(1));
      q.setZ(sn * delta(2));
    }

    return q;
  }

  static Eigen::Matrix<T,3,1> log(const Quaternion<T>& q)
  {
    // get magnitude of complex portion
    const Eigen::Matrix<T,3,1> qbar(q.x(), q.y(), q.z());
    const T qbar_mag = qbar.norm();

    // avoid numerical error with approximation
    Eigen::Matrix<T,3,1> delta;
    if (qbar_mag < T(1e-6))
    {
      if (q.w() >= T(0.0))
        delta = qbar;
      else
        delta = -qbar;
    }
    else
    {
      // Quaternions have double cover issues so wrap the axis-angle magnitude
      // to ensure the shortest rotation.
      const T delta_mag = wrapAngle(T(2.0) * atan2(qbar_mag, q.w()), T(M_PI));
      delta = delta_mag * qbar / qbar_mag;
    }

    return delta;
  }

  // q1 - q2
  static Eigen::Matrix<T,2,1> log_uvec(const Quaternion<T>& q1, const Quaternion<T>& q2)
  {
    // get unit vectors
    const Eigen::Matrix<T,3,1> e1 = q1.uvec();
    const Eigen::Matrix<T,3,1> e2 = q2.uvec();

    // avoid too small of angles
    const T e1T_e2 = e1.dot(e2);
    if (fabs(e1T_e2 - T(1.0)) < T(1e-14)) // same direction
      return Eigen::Matrix<T,2,1>(T(0.0), T(0.0));
    else if (fabs(e1T_e2 + T(1.0)) < T(1e-14)) // opposite direction
      return Eigen::Matrix<T,2,1>(T(M_PI), T(0.0));
    else
    {
      // compute axis angle difference
      const Eigen::Matrix<T,3,1> e1_x_e2 = e1.cross(e2);
      const Eigen::Matrix<T,3,1> aa = acos(e1T_e2) * e1_x_e2.normalized();

      // place error on first vector's tangent space
      return q1.proj().transpose() * aa;
    }
  }

  // derivative of quaternion exponential map
  static Eigen::Matrix<T,3,3> dexp(const Eigen::Matrix<T,3,1>& delta)
  {
    const T dmag = delta.norm();
    const Eigen::Matrix<T,3,3> delta_x = skew(delta);
    if (dmag < T(1e-6))
      return I_3x3 - T(0.5) * delta_x;
    else
    {
      const T dmag2 = dmag * dmag;
      return I_3x3 - (T(1.0) - cos(dmag)) / dmag2 * delta_x +
             (dmag - sin(dmag)) / (dmag2 * dmag) * delta_x * delta_x;
    }
  }

  T w() const { return arr(0); }
  T x() const { return arr(1); }
  T y() const { return arr(2); }
  T z() const { return arr(3); }
  void setW(const T& w) { arr(0) = w; }
  void setX(const T& x) { arr(1) = x; }
  void setY(const T& y) { arr(2) = y; }
  void setZ(const T& z) { arr(3) = z; }
  T* data() { return arr.data(); }

private:

  Eigen::Matrix<T,4,1> arr;

};

typedef Quaternion<float> Quaternionf;
typedef Quaternion<double> Quaterniond;


template<typename T = double>
class Transform
{

enum
{
  PX, PY, PZ, QW, QX, QY, QZ, T_SIZE
};

enum
{
  UX, UY, UZ, WX, WY, WZ, DT_SIZE
};

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Transform()
  {
    arr.setZero();
    arr(QW) = 1;
  }

  Transform(const Eigen::Matrix<T,3,1>& p, const Quaternion<T>& q)
  {
    setP(p);
    setQ(q);
  }

  Transform(const Eigen::Matrix<T,3,1>& p, const Eigen::Matrix<T,4,1>& q)
  {
    setP(p);
    setQ(q);
  }

  Transform(const Eigen::Matrix<T,T_SIZE,1>& t)
  {
    setT(t);
  }

  template<typename T2>
  Transform<T2> operator*(const Transform<T2>& t2)
  {
    Transform<T2> t;
    t.setP(q().rot(t2.p()) + p());
    t.setQ(t2.q() * q());
    return t;
  }

  template<typename T2>
  Transform<T2> operator+(const Eigen::Matrix<T2,DT_SIZE,1>& delta)
  {
    return *this * exp(delta);
  }

  template<typename T2>
  Eigen::Matrix<T2,DT_SIZE,1> operator-(const Transform<T2>& t1)
  {
    return log(t1.inv() * *this);
  }

  Transform<T> inv() const
  {
    return Transform<T>(-q().inv().rot(p()), q().inv());
  }

  static Transform<T> exp(const Eigen::Matrix<T,DT_SIZE,1>& delta)
  {
    Eigen::Matrix<T,3,1> u = delta.template segment<3>(UX);
    Eigen::Matrix<T,3,1> theta_vec = delta.template segment<3>(WX);
    T theta = theta_vec.norm();

    Transform<T> t;
    if (theta < T(1e-6)) // avoid numerical error with approximation
    {
      Eigen::Matrix<T,3,1> theta_vec_2 = theta_vec / T(2.0);
      t.setP(u);
      t.setQ(Eigen::Matrix<T,4,1>(1, theta_vec_2(0), theta_vec_2(1), theta_vec_2(2)));
    }
    else
    {
      T theta2 = theta * theta;
      T theta3 = theta * theta2;
      Eigen::Matrix<T,3,3> theta_skew = skew(theta_vec);
      Eigen::Matrix<T,3,3> theta_skew2 = theta_skew * theta_skew;
      Eigen::Matrix<T,3,3> V = I_3x3.cast<T>() + (T(1.0) - cos(theta))/theta2 *
                               theta_skew + (theta - sin(theta))/theta3 * theta_skew2;
      t.setP(V * u);
      t.setQ(Quaternion<T>::exp(theta_vec));
    }

    return t;
  }

  static Eigen::Matrix<T,DT_SIZE,1> log(const Transform<T>& t)
  {
    Eigen::Matrix<T,3,1> theta_vec = Quaternion<T>::log(t.q());
    T theta = theta_vec.norm();

    // avoid numerical error with approximation
    Eigen::Matrix<T,DT_SIZE,1> delta;
    if (theta < T(1e-6))
    {
      delta.template segment<3>(PX) = t.p();
    }
    else
    {
      T theta2 = theta * theta;
      Eigen::Matrix<T,3,3> theta_skew = skew(theta_vec);
      Eigen::Matrix<T,3,3> theta_skew2 = theta_skew * theta_skew;
      Eigen::Matrix<T,3,3> Vinv = I_3x3.cast<T>() - T(0.5) * theta_skew + (T(1.0)/theta2) * (T(1.0) -
                                  theta*sin(theta)/(T(2.0)*(T(1.0) - cos(theta)))) * theta_skew2;
      delta.template segment<3>(PX) = Vinv * t.p();
    }
    delta.template segment<3>(WX) = theta_vec;

    return delta;
  }

  void setT(const Eigen::Matrix<T,T_SIZE,1> t) { arr = t; }
  void setP(const Eigen::Matrix<T,3,1> p) { arr.template segment<3>(PX) = p; }
  void setQ(const Quaternion<T> q) { arr.template segment<4>(QW) = q.toEigen(); }
  void setQ(const Eigen::Matrix<T,4,1> q) { arr.template segment<4>(QW) = q; }
  void setPX(const T& x) { arr(PX) = x; }
  void setPY(const T& y) { arr(PY) = y; }
  void setPZ(const T& z) { arr(PZ) = z; }
  void setQW(const T& w) { arr(QW) = w; }
  void setQX(const T& x) { arr(QX) = x; }
  void setQY(const T& y) { arr(QY) = y; }
  void setQZ(const T& z) { arr(QZ) = z; }

  Eigen::Matrix<T,3,1> p() const { return arr.template segment<3>(PX); }
  Quaternion<T> q() const { return Quaternion<T>(arr(QW), arr(QX), arr(QY), arr(QZ)); }
  T* data() const { return arr.data(); }

private:

  Eigen::Matrix<T,T_SIZE,1> arr;

};

typedef Transform<float> Transformf;
typedef Transform<double> Transformd;


} // namespace common

#endif // COMMON_H
