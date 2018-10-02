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

  Quaternion()
  {
    w = T(1.0);
    x = T(0.0);
    y = T(0.0);
    z = T(0.0);
  }

  Quaternion(const T& _w, const T& _x, const T& _y, const T& _z)
  {
    w = _w;
    x = _x;
    y = _y;
    z = _z;
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
    w = cr*cp*cy + sr*sp*sy;
    x = sr*cp*cy - cr*sp*sy;
    y = cr*sp*cy + sr*cp*sy;
    z = cr*cp*sy - sr*sp*cy;
  }

  Quaternion(const Eigen::Matrix<T,4,1>& v)
  {
    w = v(0);
    x = v(1);
    y = v(2);
    z = v(3);
  }

  // create quaternion from a unit vector
  Quaternion(Eigen::Matrix<T,3,1>& fz)
  {
    // convert to axis-angle representation
    fz.normalize(); // enforce unit length
    const T theta = acos(fz.dot(e1.cast<T>()));

    if (theta < T(1e-6))
    {
      w = T(1.0);
      x = T(0.0);
      y = T(0.0);
      z = T(0.0);
    }
    else
    {
      const Eigen::Matrix<T,3,1> iaa = (fz.cross(e1.cast<T>())).normalized();

      // get complex portion of quaternion
      const T theta_2 = theta / T(2.0);
      const Eigen::Matrix<T,3,1> qv = iaa * sin(theta_2);

      w = cos(theta_2);
      x = qv(0);
      y = qv(1);
      z = qv(2);
    }
  }

  void normalize()
  {
    const T m = mag();
    w /= m;
    x /= m;
    y /= m;
    z /= m;
  }

  // initialize random unit quaternion
  Quaternion(std::normal_distribution<T>& dist, std::default_random_engine& rng)
  {
    Quaternion q(dist(rng), dist(rng), dist(rng), dist(rng));
    q.normalize();
    q.w = (q.w > 0) ? q.w : -q.w;
    w = q.w;
    x = q.x;
    y = q.y;
    z = q.z;
  }

  // overload addition operator as boxplus for a quaternion and a 3-vector
  Quaternion operator+(const Eigen::Matrix<T,3,1>& delta) const
  {
    return *this * exp(delta);
  }

  void operator+=(const Eigen::Matrix<T,3,1>& delta)
  {
    *this = *this + delta;
  }

  // overload minus operator as boxminus for two quaternions
  Eigen::Matrix<T,3,1> operator-(const Quaternion &q2) const
  {
    return log(q2.inv() * *this);
  }

  Quaternion operator*(const Quaternion &q2) const
  {
    const T qw = w*q2.w - x*q2.x - y*q2.y - z*q2.z;
    const T qx = w*q2.x + x*q2.w + y*q2.z - z*q2.y;
    const T qy = w*q2.y - x*q2.z + y*q2.w + z*q2.x;
    const T qz = w*q2.z + x*q2.y - y*q2.x + z*q2.w;
    return Quaternion(qw, qx, qy, qz);
  }

  void operator*=(const Quaternion &q)
  {
    *this = *this * q;
  }

  friend std::ostream& operator<<(std::ostream &os, const Quaternion &q)
  {
    os << q.w << "\n" << q.x << "\n" << q.y << "\n" << q.z << "\n";
    return os;
  }

  void scale(const T& s)
  {
    w *= s;
    x *= s;
    y *= s;
    z *= s;
  }

  T mag() const
  {
    return sqrt(w*w + x*x + y*y + z*z);
  }

  T roll() const
  {
    return atan2(T(2.0) * (w*x + y*z), T(1.0) - T(2.0) * (x*x + y*y));
  }

  T pitch() const
  {
    const T val = T(2.0) * (w*y - x*z);

    // hold at 90 degrees if invalid
    if (fabs(val) > T(1.0))
      return copysign(T(1.0), val) * T(M_PI) / T(2.0);
    else
      return asin(val);
  }

  T yaw() const
  {
    return atan2(T(2.0) * (w*z + x*y), T(1.0) - T(2.0) * (y*y + z*z));
  }

  Eigen::Matrix<T,3,1> euler() const
  {
    return Eigen::Matrix<T,3,1>(roll(),pitch(),yaw());
  }

  void fromEigen(const Eigen::Matrix<T,4,1>& q)
  {
    w = q(0);
    x = q(1);
    y = q(2);
    z = q(3);
  }

  Eigen::Matrix<T,4,1> toEigen() const
  {
    return Eigen::Matrix<T,4,1>(w,x,y,z);
  }

  Eigen::Matrix<T,3,1> bar() const
  {
    return Eigen::Matrix<T,3,1>(x,y,z);
  }

  Quaternion inv() const
  {
    return Quaternion<T>(w, -x, -y, -z);
  }

  Eigen::Matrix<T,3,3> R() const
  {
    // Pre-calculations
    const T ww = w * w;
    const T wx = w * x;
    const T wy = w * y;
    const T wz = w * z;
    const T xy = x * y;
    const T xz = x * z;
    const T yz = y * z;

    // Output
    Eigen::Matrix<T,3,3> R;
    R(0,0) = T(2.0) * (ww + x * x) - T(1.0);
    R(0,1) = T(2.0) * (wz + xy);
    R(0,2) = T(2.0) * (xz - wy);
    R(1,0) = T(2.0) * (xy - wz);
    R(1,1) = T(2.0) * (ww + y * y) - T(1.0);
    R(1,2) = T(2.0) * (wx + yz);
    R(2,0) = T(2.0) * (wy + xz);
    R(2,1) = T(2.0) * (yz - wx);
    R(2,2) = T(2.0) * (ww + z * z) - T(1.0);
    return R;
  }

  Eigen::Matrix<T,3,1> rotSlow(const Eigen::Matrix<T,3,1>& v) const
  {
    const Quaternion qv(T(0.0), v(0), v(1), v(2));
    const Quaternion qv_new = inv() * qv * *this;
    return Eigen::Matrix<T,3,1>(qv_new.x, qv_new.y, qv_new.z);
  }

  Eigen::Matrix<T,3,1> rot(const Eigen::Matrix<T,3,1>& v) const
  {
    const Eigen::Matrix<T,3,1> t = T(2.0) * v.cross(bar());
    return v + w * t + t.cross(bar());
  }

  Eigen::Matrix<T,3,1> uvec() const
  {
    return rot(e3.cast<T>());
  }

  Eigen::Matrix<T,3,2> proj() const
  {
    return R() * I_2x3.cast<T>().transpose();
  }

  static Quaternion exp(const Eigen::Matrix<T,3,1>& delta)
  {
    const T delta_norm = delta.norm();

    Quaternion<T> q;
    if (delta_norm < T(1e-6)) // avoid numerical error with approximation
    {
      q.w = T(1.0);
      q.x = delta(0) / T(2.0);
      q.y = delta(1) / T(2.0);
      q.z = delta(2) / T(2.0);
    }
    else
    {
      const T delta_norm_2 = delta_norm / T(2.0);
      const T sn = sin(delta_norm_2) / delta_norm;
      q.w = cos(delta_norm_2);
      q.x = sn * delta(0);
      q.y = sn * delta(1);
      q.z = sn * delta(2);
    }

    return q;
  }

  static Eigen::Matrix<T,3,1> log(const Quaternion& q)
  {
    // get magnitude of complex portion
    const Eigen::Matrix<T,3,1> qbar(q.x, q.y, q.z);
    const T qbar_mag = qbar.norm();

    // avoid numerical error with approximation
    Eigen::Matrix<T,3,1> delta;
    if (qbar_mag < T(1e-6))
    {
      if (q.w >= T(0.0))
        delta = qbar;
      else
        delta = -qbar;
    }
    else
    {
      // Quaternions have double cover issues so wrap the axis-angle magnitude
      // to ensure the shortest rotation.
      const T delta_mag = wrapAngle(T(2.0) * atan2(qbar_mag, q.w), T(M_PI));
      delta = delta_mag * qbar / qbar_mag;
    }

    return delta;
  }

  // q1 - q2
  static Eigen::Matrix<T,2,1> log_uvec(const Quaternion& q1, const Quaternion& q2)
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

  T getW() const { return w; }
  T getX() const { return x; }
  T getY() const { return y; }
  T getZ() const { return z; }

private:

  T w;
  T x;
  T y;
  T z;

};


} // namespace common

#endif // COMMON_H
