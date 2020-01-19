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


// Approximate constants near Earth's surface
static const double gravity = 9.80665; // (m/s^2) Gravitational acceleration
static const double R_earth = 6371008.8; // (m) Earth's mean radius
static const double B0 = 31200.0; // (nT) mean magnetic field at Earth's equator
static const double P_sea = 101325.0; // (Pa) Standard atmospheric pressure at sea level
static const double T_sea = 288.15; // (K) Standard temperature at sea level
static const double T_lapse = 0.0065; // (K/m) Temperature lapse rate
static const double R_gas = 8.31447; // (J/(mol*K) Universal gas constant
static const double M_air = 0.0289644; // (kg/mol) Molar mass of dry air
static const double MNP_lat = 80.37; // (deg) magnetic north pole latitude as of 2015
static const double MNP_lon = -72.62; // (deg) magnetic north pole longitude as of 2015


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


// Atmospheric pressure as a function of altitude above sea level
template<typename T>
T airPres(const T &alt)
{
  return P_sea * pow(1.0 - T_lapse * alt / T_sea, gravity * M_air / (R_gas * T_lapse));
}


// Air density as a function of altitude (meters) above sea level and air temperature (degrees Fahrenheit)
template<typename T>
T airDense(const T &alt, const T &temp)
{
  T tempK = 5.0 / 9.0 * (temp - 32.0) + 273.15; // convert to Kelvin
  return airPres(alt) * M_air / (R_gas * tempK);
}


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
T round2dec(const T &number, const int &decimal_place)
{
  T val = pow(T(10),decimal_place);
  return round(number * val) / val;
}


// skew symmetric matrix from vector
template<typename T>
Eigen::Matrix<typename T::Scalar,3,3> skew(const Eigen::MatrixBase<T>& vec)
{
  Eigen::Matrix<typename T::Scalar,3,3> A;
  typename T::Scalar zr(0);
  A <<      zr, -vec(2),  vec(1),
        vec(2),      zr, -vec(0),
       -vec(1),  vec(0),      zr;
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
  if (theta > 1e-6)
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
      printf("Eigen Matrix Size does not match parameter size for \"%s\" in %s\n", key.c_str(), filename.c_str());
      return false;
    }
  }
  else
  {
    printf("Unable to load \"%s\" from %s\n", key.c_str(), filename.c_str());
    return false;
  }
}


// Load array from a .yaml file into a diagonal Eigen-type matrix
template <typename T>
bool get_yaml_eigen_diag(const std::string key, const std::string filename, Eigen::MatrixBase<T>& val)
{
  Eigen::Matrix<typename T::Scalar, T::RowsAtCompileTime, 1> val2;
  get_yaml_eigen(key, filename, val2);
  val = val2.asDiagonal();
  return true;
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


template<typename T, int S>
Eigen::Matrix<T,S,1> saturateVector(const T& max, Eigen::Matrix<T,S,1>& vec)
{
  T vec_mag = vec.norm();
  if (vec_mag <= max)
    return vec;
  else
    return vec / vec_mag * max;
}


// Random normal matrix generation
template<typename T>
void randomNormal(Eigen::DenseBase<T>& matrix,
                  std::normal_distribution<typename T::Scalar>& dist,
                  std::default_random_engine& rng)
{
  for (int i = 0; i < matrix.rows(); ++i)
    for (int j = 0; j < matrix.cols(); ++j)
      matrix(i,j) = dist(rng);
}


// Perspective projection into an image
template<typename T>
void projToImg(Eigen::Matrix<T,2,1>& pix, const Eigen::Matrix<T,3,1> &lc, const Eigen::Matrix<T,3,3> &K)
{
  pix = K.topRows(2) * (lc / lc(2));
}


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


// Sign function
template<typename T>
int sign(T val)
{
    return (T(0) <= val) - (val < T(0));
}


// Generic PID controller class
template<typename T = double>
struct PID
{
  PID()
  {
    kp = T(0);
    ki = T(0);
    kd = T(0);
    max = T(0);
    integrator = T(0);
    differentiator = T(0);
    prev_x = T(0);
    tau = T(0);
  }

  void init(T _kp, T _ki, T _kd, T _max, T _min, T _tau)
  {
    kp = _kp;
    ki = _ki;
    kd = _kd;
    max = _max;
    min = _min;
    tau = _tau;
  }

  T run(T dt, T x, T x_c, bool update_integrator)
  {
    T xdot;
    if (dt > 1e-4)
    {
      // calculate D term (use dirty derivative if we don't have access to a measurement of the derivative)
      // The dirty derivative is a sort of low-pass filtered version of the derivative.
      //// (Include reference to Dr. Beard's notes here)
      differentiator = (T(2) * tau - dt) / (T(2) * tau + dt) * differentiator
          + T(2) / (T(2) * tau + dt) * (x - prev_x);
      xdot = differentiator;
    }
    else
    {
      xdot = T(0);
    }
    prev_x = x;

    return run(dt, x, x_c, update_integrator, xdot);
  }

  T run(T dt, T x, T x_c, bool update_integrator, T xdot)
  {
    // Calculate Error
    T error = x_c - x;

    // Initialize Terms
    T p_term = error * kp;
    T i_term = T(0);
    T d_term = T(0);

    // If there is a derivative term
    if (kd > T(0))
    {
      d_term = kd * xdot;
    }

    //If there is an integrator term and we are updating integrators
    if ((ki > T(0)) && update_integrator)
    {
      // integrate
      integrator += error * dt;
      // calculate I term
      i_term = ki * integrator;
    }

    // sum three terms
    T u = p_term - d_term + i_term;

    // Integrator anti-windup
    T u_sat = (u > max) ? max : (u < min) ? min : u;
    if (u != u_sat && std::abs(i_term) > std::abs(u - p_term + d_term) && ki > T(0))
      integrator = (u_sat - p_term + d_term)/ki;

    // Set output
    return u_sat;
  }

  T kp, ki, kd;
  T max, min;
  T integrator, differentiator;
  T prev_x;
  T tau;
};


// Solution to general cubic equation
// https://en.wikipedia.org/wiki/Cubic_equation
inline void solveCubic(const double& _a, const double& _b, const double& _c, const double& _d, std::vector<double>& real, std::vector<double>& imag)
{
  // Make coefficients complex
  std::complex<double> a(_a, 0.0);
  std::complex<double> b(_b, 0.0);
  std::complex<double> c(_c, 0.0);
  std::complex<double> d(_d, 0.0);

  // Constant
  static std::complex<double> xi = (-1.0 + sqrt(std::complex<double>(-3.0)))/2.0;

  // Compute 3 possible solutions
  std::complex<double> Delta_0 = b*b - 3.0*a*c;
  std::complex<double> Delta_1 = 2.0*b*b*b - 9.0*a*b*c + 27.0*a*a*d;

  std::complex<double> dd = sqrt(pow(Delta_1,2.0) - 4.0*pow(Delta_0,3.0));
  std::complex<double> C = pow((Delta_1 + dd)/2.0, 1.0/3.0);
  if (std::abs(C.real()) < 1e-6)
    C = pow((Delta_1 - dd)/2.0, 1.0/3.0);
  
  std::vector<std::complex<double> > solutions;
  for (int k = 0; k < 3; ++k)
    solutions.push_back(-1.0/3.0*(b + pow(xi,double(k))*C + Delta_0/(pow(xi,double(k))*C)));
  
  // Sort solutions in ascending order of real parts
  std::sort(solutions.begin(), solutions.end(),
            [](const std::complex<double>& a, const std::complex<double>& b)
            { return a.real() < b.real(); });

  // Pack output
  real.clear();
  imag.clear();
  for (const auto& s : solutions)
  {
    real.push_back(s.real());
    imag.push_back(s.imag());
  }
}


} // namespace common

#endif // COMMON_H
