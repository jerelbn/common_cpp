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

typedef Eigen::Matrix<double, 1, 1> Vector1d;
typedef Eigen::Matrix<double, 5, 1> Vector5d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 7, 1> Vector7d;
typedef Eigen::Matrix<double, 8, 1> Vector8d;
typedef Eigen::Matrix<double, 9, 1> Vector9d;

typedef Eigen::Matrix<double, 1, 1> Matrix1d;
typedef Eigen::Matrix<double, 5, 5> Matrix5d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 7, 7> Matrix7d;
typedef Eigen::Matrix<double, 8, 8> Matrix8d;
typedef Eigen::Matrix<double, 9, 9> Matrix9d;

typedef Eigen::Matrix<float, 1, 1> Vector1f;
typedef Eigen::Matrix<float, 5, 1> Vector5f;
typedef Eigen::Matrix<float, 6, 1> Vector6f;
typedef Eigen::Matrix<float, 7, 1> Vector7f;
typedef Eigen::Matrix<float, 8, 1> Vector8f;
typedef Eigen::Matrix<float, 9, 1> Vector9f;

typedef Eigen::Matrix<float, 1, 1> Matrix1f;
typedef Eigen::Matrix<float, 5, 5> Matrix5f;
typedef Eigen::Matrix<float, 6, 6> Matrix6f;
typedef Eigen::Matrix<float, 7, 7> Matrix7f;
typedef Eigen::Matrix<float, 8, 8> Matrix8f;
typedef Eigen::Matrix<float, 9, 9> Matrix9f;

namespace common
{


// Approximate constants near Earth's surface
static constexpr double gravity = 9.80665;   // (m/s^2) Gravitational acceleration
static constexpr double R_earth = 6371008.8; // (m) Earth's mean radius
static constexpr double B0      = 31200.0;   // (nT) mean magnetic field at Earth's equator
static constexpr double P_sea   = 101325.0;  // (Pa) Standard atmospheric pressure at sea level
static constexpr double T_sea   = 288.15;    // (K) Standard temperature at sea level
static constexpr double T_lapse = 0.0065;    // (K/m) Temperature lapse rate
static constexpr double R_gas   = 8.31447;   // (J/(mol*K) Universal gas constant
static constexpr double M_air   = 0.0289644; // (kg/mol) Molar mass of dry air
static constexpr double MNP_lat = 80.37;     // (deg) magnetic north pole latitude as of 2015
static constexpr double MNP_lon = -72.62;    // (deg) magnetic north pole longitude as of 2015

// Constant vectors and matrices
static const Eigen::Vector3d e1(1, 0, 0);
static const Eigen::Vector3d e2(0, 1, 0);
static const Eigen::Vector3d e3(0, 0, 1);
static const Eigen::Matrix3d I_3x3 = Eigen::Matrix3d::Identity();
static const Eigen::Matrix2d I_2x2 = Eigen::Matrix2d::Identity();
static const Eigen::Matrix<double,2,3> I_2x3 = (Eigen::Matrix<double,2,3>() << 1, 0, 0, 0, 1, 0).finished();
static const Eigen::Matrix<double,3,2> I_3x2 = (Eigen::Matrix<double,3,2>() << 1, 0, 0, 1, 0, 0).finished();


/** @brief Computes atmospheric pressure as a function of altitude.
 @param alt altitude above sea level in meters
 */
template<typename T>
T airPres(const T &alt)
{
  return P_sea * pow(1.0 - T_lapse * alt / T_sea, gravity * M_air / (R_gas * T_lapse));
}


/** @brief Computes air density as a function of altitude and temperature.
 @param alt altitude above sea level in meters
 @param temp air temperature in degrees Fahrenheit
 */
template<typename T>
T airDense(const T &alt, const T &temp)
{
  T tempK = 5.0 / 9.0 * (temp - 32.0) + 273.15; // convert to Kelvin
  return airPres(alt) * M_air / (R_gas * tempK);
}


/** @brief Wraps an angle to some bound.
 @param angle angle to be wrapped in radians
 @param bound boundary of wrapping (usually 2*pi or pi)
 */
template<typename T>
T wrapAngle(const T &angle, const T &bound)
{
  if (angle > bound)
    return angle - T(2.0) * T(M_PI);
  if (angle < bound - T(2.0) * T(M_PI))
    return angle + T(2.0) * T(M_PI);
  return angle;
}


/** @brief Returns the sign of input value.
 @param val a real number
 */
template<typename T>
int sign(T val)
{
    return (T(0) <= val) - (val < T(0));
}


/** Round to a decimal place.
 @param number number to round
 @param decimal_place integer value of decimal place to round to
 */
template<typename T>
T roundDecimal(const T &number, const int &decimal_place)
{
  T val = pow(T(10),decimal_place);
  return round(number * val) / val;
}


/** Create skew symmetric matrix from a vector.
 @param vec vector from which to create matrix
 */
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


/** @brief Creates a vector from a skew symmetric matrix.
 @param mat skew symmetric matrix from which to create vector
 */
template<typename T>
Eigen::Matrix<T,3,1> vex(const Eigen::Matrix<T,3,3>& mat)
{
  return (Eigen::Matrix<T,3,1>() << mat(2,1), mat(0,2), mat(1,0)).finished();
}


/** @brief Angular difference between two vectors.
 @param v1 first vector
 @param v2 second vector
 */
template<typename T>
T angleBetweenVectors(const Eigen::Matrix<T,3,1>& v1, const Eigen::Matrix<T,3,1>& v2)
{
  T val = (v1.transpose() * v2)(0) / (v1.norm() * v2.norm());
  if (val > 1)
    return 0;
  else if (val < -1)
    return M_PI;
  else
    return acos(val);
}


/** @brief Saturates a scalar value.
 @param val number to be saturated
 @param max maximum allowed value
 @param min minimum allowed value
 */
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


/** @brief Create unit vector from camera image pixel coordinates and intrinsic camera matrix.
 @param uvec output unit vector in camera coordinates
 @param pix pixel position in image
 @param K camera intrinsic matrix
 */
template<typename T>
void unitVectorFromPixelPosition(Eigen::Matrix<T,3,1>& uvec, const Eigen::Matrix<T,2,1> &pix, const Eigen::Matrix<T,3,3> &K)
{
  uvec(0) = (pix(0) - K(0,2))/K(0,0);
  uvec(1) = (pix(1) - K(1,2))/K(1,1);
  uvec(2) = T(1.0);
  uvec.normalize();
}


/** @brief Perspective projection of vector in camera coordinates into image coordinates.
 @param pix output of pixel position in image
 @param lc vector in camera coordinates
 @param K camera intrinsic matrix
 */
template<typename T>
void projectToImage(Eigen::Matrix<T,2,1>& pix, const Eigen::Matrix<T,3,1> &lc, const Eigen::Matrix<T,3,3> &K)
{
  pix = K.topRows(2) * (lc / lc(2));
}


/** @brief Computes unit vector from pixel coordinates.
 @param dir output of unit vector
 @param pix pixel position in image coordinates
 @param K_inv inverse of camera intrinsic matrix
 */
template<typename T>
void dirFromPix(Eigen::Matrix<T,3,1> &dir, const Eigen::Matrix<T,2,1> &pix, const Eigen::Matrix<T,3,3> &K_inv)
{
  dir = K_inv * Eigen::Matrix<T,3,1>(pix(0), pix(1), 1);
  dir /= dir.norm();
}


/** @brief Loads scalar parameters from a .yaml file.
 @note Original author: James Jackson
 @param key string of parameter name
 @param filename name of .yaml file to load parameters from
 @param val output value of loaded parameter
 @param print_error error printing boolean - set to false to suppress error printing
 */
template <typename T>
bool getYamlNode(const std::string& key, const std::string& filename, T& val, const bool& print_error = true)
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


/** @brief Loads array from a .yaml file into an Eigen-type matrix or vector.
 @note Original author: James Jackson
 @param key string of parameter name
 @param filename name of .yaml file to load parameters from
 @param val output matrix or vector of loaded parameters
 */
template <typename T>
bool getYamlEigen(const std::string key, const std::string filename, Eigen::MatrixBase<T>& val)
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


/** @brief Loads array from a .yaml file into an Eigen-type diagonal matrix.
 @note Original author: James Jackson
 @param key string of parameter name
 @param filename name of .yaml file to load parameters from
 @param val output matrix of loaded parameters
 */
template <typename T>
bool getYamlEigenDiag(const std::string key, const std::string filename, Eigen::MatrixBase<T>& val)
{
  Eigen::Matrix<typename T::Scalar, T::RowsAtCompileTime, 1> val2;
  getYamlEigen(key, filename, val2);
  val = val2.asDiagonal();
  return true;
}


/** @brief Loads array from a binary file and returns the pointer to the beginning of the array.
 @param filename name of binary file to load parameters from
 @param array_size size of array stored in binary file
 */
template<typename T>
T* loadBinary(const std::string& filename, long& array_size)
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


/** @brief Loads binary file directly into an Eigen-type matrix.
 @param filename name of binary file to load parameters from
 @param matrix_rows number of rows in output matrix
 */
template<typename T>
Eigen::Matrix<T,-1,-1> loadBinaryToMatrix(const std::string& filename, const int& matrix_rows)
{
  long array_size;
  T* ptr = loadBinary<T>(filename, array_size);
  Eigen::Matrix<T,-1,-1> m(matrix_rows, array_size/matrix_rows);
  copy_ptr_to_eigen(ptr, m);
  delete[] ptr;
  return m;
}


/** @brief Loads binary file directly into a standard vector
 @param filename name of binary file to load parameters from
 */
template<typename T>
std::vector<T> loadBinaryToVector(const std::string& filename)
{
  long array_size;
  T* ptr = loadBinary<T>(filename, array_size);
  std::vector<T> v(array_size, T(0));
  for (int i = 0; i < v.size(); ++i)
    v[i] = *(ptr+i);
  delete[] ptr;
  return v;
}


/** @brief A generic PID controller.
 * @param kp proportional gain
 * @param ki integral gain
 * @param kd derivative gain
 * @param max maximum limit on output
 * @param min minimum limit on output
 * @param integrator storage for integrated error
 * @param differentiator storage for derivative in numerical derivative
 * @param prev_x previous state
 * @param tau time constant for numerical derivative
 */
template<typename T = double>
struct PID
{
  PID()
  {
    kp = T(0);
    ki = T(0);
    kd = T(0);
    umax = T(0);
    umin = T(0);
    integrator = T(0);
    differentiator = T(0);
    prev_x = T(0);
    tau = T(0);
  }

  void init(T _kp, T _ki, T _kd, T _umax, T _umin, T _emax, T _emin, T _tau)
  {
    kp = _kp;
    ki = _ki;
    kd = _kd;
    umax = _umax;
    umin = _umin;
    emax = _emax;
    emin = _emin;
    tau = _tau;
  }

  T run(T dt, T x, T x_c, bool update_integrator)
  {
    T xdot;
    if (dt > T(1e-8))
    {
      // Low pass filtered derivative
      T tau_x_2 = T(2) * tau;
      T tau_x_2_plus_dt = tau_x_2 + dt;
      differentiator = (tau_x_2 - dt) / tau_x_2_plus_dt * differentiator + T(2) / tau_x_2_plus_dt * (x - prev_x);
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
    // Calculate error and initialize components
    T error = x_c - x;
    T error_sat = saturate(error, emax, emin);
    T p_term = kp * error_sat;
    T i_term = T(0);
    T d_term = T(0);

    // Compute ID components and sum them for total output
    if (kd > T(0))
      d_term = kd * xdot;

    if (ki > T(0) && update_integrator)
    {
      integrator += error * dt;
      i_term = ki * integrator;
    }

    T u = p_term - d_term + i_term;

    // Integrator anti-windup
    T u_sat = saturate(u, umax, umin);
    if (ki > T(0) && (u != u_sat || error != error_sat))
      integrator = dt/ki*(u_sat - u);

    return u_sat;
  }

  T kp, ki, kd;
  T umax, umin;
  T emax, emin;
  T integrator, differentiator;
  T prev_x;
  T tau;
};


/** @brief Computes solution to the general cubic equation (see https://en.wikipedia.org/wiki/Cubic_equation).
 * @param _a first coefficient of cubic equation
 * @param _b second coefficient of cubic equation
 * @param _c third coefficient of cubic equation
 * @param _d fourth coefficient of cubic equation
 * @param real output of real parts of solution
 * @param imag output of complex parts of solution
 */
template<typename T = double>
inline void solveCubic(const T& _a, const T& _b, const T& _c, const T& _d, std::vector<T>& real, std::vector<T>& imag)
{
  // Make coefficients complex
  std::complex<T> a(_a, T(0.0));
  std::complex<T> b(_b, T(0.0));
  std::complex<T> c(_c, T(0.0));
  std::complex<T> d(_d, T(0.0));

  // Constant
  static std::complex<T> xi = (-T(1.0) + sqrt(std::complex<T>(-T(3.0))))/T(2.0);

  // Compute 3 possible solutions
  std::complex<T> Delta_0 = b*b - T(3.0)*a*c;
  std::complex<T> Delta_1 = T(2.0)*b*b*b - 9.0*a*b*c + T(27.0)*a*a*d;

  std::complex<T> dd = sqrt(pow(Delta_1,T(2.0)) - T(4.0)*pow(Delta_0,T(3.0)));
  std::complex<T> C = pow((Delta_1 + dd)/T(2.0), T(1.0)/T(3.0));
  if (std::abs(C.real()) < T(1e-8))
    C = pow((Delta_1 - dd)/T(2.0), T(1.0)/T(3.0));
  
  std::vector<std::complex<T> > solutions;
  for (int k = 0; k < 3; ++k)
    solutions.push_back(-T(1.0)/T(3.0)*(b + pow(xi,T(k))*C + Delta_0/(pow(xi,T(k))*C)));
  
  // Sort solutions in ascending order of real parts
  std::sort(solutions.begin(), solutions.end(),
            [](const std::complex<T>& a, const std::complex<T>& b)
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


/** @brief Rotation from camera body to camera coordinates (Front-Right-Down --> Right-Down-Front)
 */
static const Eigen::Matrix3d R_cb2c = (Eigen::Matrix3d() << 0, 1, 0, 0, 0, 1, 1, 0, 0).finished();


/** @brief Returns rotation matrix for positive coordinate frame rotation about the x-axis.
 @param roll rotation about x-axis
 */
template<typename T>
Eigen::Matrix<T,3,3> Rx(const T& roll)
{
  return (Eigen::Matrix<T,3,3>() << 1,          0,         0,
                                    0,  cos(roll), sin(roll),
                                    0, -sin(roll), cos(roll)).finished();
}

/** @brief Returns rotation matrix for positive coordinate frame rotation about the y-axis.
 @param pitch rotation about y-axis
 */
template<typename T>
Eigen::Matrix<T,3,3> Ry(const T& pitch)
{
  return (Eigen::Matrix<T,3,3>() << cos(pitch), 0, -sin(pitch),
                                             0, 1,           0,
                                    sin(pitch), 0,  cos(pitch)).finished();
}


/** @brief Returns rotation matrix for positive coordinate frame rotation about the z-axis.
 @param yaw rotation about z-axis
 */
template<typename T>
Eigen::Matrix<T,3,3> Rz(const T& yaw)
{
  return (Eigen::Matrix<T,3,3>() <<  cos(yaw), sin(yaw), 0,
                                    -sin(yaw), cos(yaw), 0,
                                            0,        0, 1).finished();
}


/** @brief Returns rotation matrix for positive coordinate frame rotation about the x-y-z axes (Euler 3-2-1 or Z-Y-X).
 @param roll rotation about x-axis
 @param pitch rotation about y-axis
 @param yaw rotation about z-axis
 */
template<typename T>
Eigen::Matrix<T,3,3> Rzyx(const T& roll, const T& pitch, const T& yaw)
{
  return Rx(roll) * Ry(pitch) * Rz(yaw);
}


/** @brief Returns rotation matrix for positive coordinate frame rotation about the x-y-z axes (Euler 3-1-2 or Z-X-Y).
 @param roll rotation about x-axis
 @param pitch rotation about y-axis
 @param yaw rotation about z-axis
 */
template<typename T>
Eigen::Matrix<T,3,3> Rzxy(const T& roll, const T& pitch, const T& yaw)
{
  return Ry(pitch) * Rx(roll) * Rz(yaw);
}


/** @brief Computes matrix exponential for a given skew symmetric matrix.
 @param deltax skew symmetric matrix input
 */
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


/** @brief Computes the logarithmic map of a rotation matrix.
 @param R rotation matrix input
 */
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


/** @brief Returns rotation matrix rotating Euler angular rates to body angular rates.
 @note This supports 3-2-1 and 3-1-2 Euler angle orders.
 @param roll rotation about body x-axis
 @param pitch rotation about body y-axis
 @param order euler angle rotation order
 */
template<typename T>
Eigen::Matrix<T,3,3> R_euler_rate_to_body_rate(const T& roll, const T& pitch, const int& order)
{
  // Pre-calcs
  double sr = sin(roll);
  double cr = cos(roll);
  double sp = sin(pitch);
  double cp = cos(pitch);

  // Populate matrix for different orders
  Eigen::Matrix<T,3,3> R = Eigen::Matrix<T,3,3>::Identity();
  if (order == 321)
  {
    R(0,2) = -sp;
    R(1,1) = cr;
    R(1,2) = sr * cp;
    R(2,1) = -sr;
    R(2,2) = cr * cp;
  }
  else if (order == 312)
  {
    R(0,0) = cp;
    R(0,2) = -cr * sp;
    R(1,2) = sr;
    R(2,0) = sp;
    R(2,2) = cr * cp;
  }
  else
  {
    std::stringstream ss;
    ss << "Invalid Euler angle order in " << __FILE__ \
        << " line " << __LINE__ << ": " << "Choose 321 or 312." << std::endl;
    throw std::runtime_error(ss.str());
  }

  return R;
}


/** @brief 4th order integrator for arbitrary size arrays
 @param f function containing equations of motion
 @param dt time step of integration
 @param x state of XSIZE
 @param u input of USIZE
 @param dx change in state of XSIZE
 */
template<typename T, int XSIZE, int USIZE>
void rk4(std::function<void(const T[XSIZE], const T[USIZE], T[XSIZE])> f,
                            const T& dt,
                            const T x[XSIZE],
                            const T u[USIZE],
                            T dx[XSIZE])
{
  int i;
  T k1[XSIZE], k2[XSIZE], k3[XSIZE], k4[XSIZE];
  T x1[XSIZE], x2[XSIZE], x3[XSIZE];
  f(x,  u, k1);
  for (i = 0; i < XSIZE; ++i)
    x1[i] = x[i] + k1[i] * dt / T(2.0);
  f(x1, u, k2);
  for (i = 0; i < XSIZE; ++i)
    x2[i] = x[i] + k2[i] * dt / T(2.0);
  f(x2, u, k3);
  for (i = 0; i < XSIZE; ++i)
    x3[i] = x[i] + k3[i] * dt;
  f(x3, u, k4);
  for (i = 0; i < XSIZE; ++i) {
    dx[i] = (k1[i] + T(2.0) * k2[i] + T(2.0) * k3[i] + k4[i]) * dt / T(6.0);
  }
}


/** @brief 4th order integrator for arbitrary size vectors of Eigen::Matrix type
 @param f function containing equations of motion
 @param dt time step of integration
 @param x state of XSIZE
 @param u input of USIZE
 @param dx change in state of XSIZE
 */
template<typename T, int XSIZE, int USIZE>
void rk4(std::function<void(const Eigen::Matrix<T,XSIZE,1>&, const Eigen::Matrix<T,USIZE,1>&, Eigen::Matrix<T,XSIZE,1>&)> f,
                            const T& dt,
                            const Eigen::Matrix<T,XSIZE,1>& x,
                            const Eigen::Matrix<T,USIZE,1>& u,
                            Eigen::Matrix<T,XSIZE,1>& dx)
{
  Eigen::Matrix<T,XSIZE,1> k1, k2, k3, k4;
  f(x, u, k1);
  f(x + k1 * dt / T(2.0), u, k2);
  f(x + k2 * dt / T(2.0), u, k3);
  f(x + k3 * dt, u, k4);
  dx = (k1 + T(2.0) * k2 + T(2.0) * k3 + k4) * dt / T(6.0);
}


/** @brief A templated circular buffer
 @param buf Data buffer
 @param max_size Size of array when full
 @param head Index to beginning of the array
 @param tail Index to the end of the array
 @param full Indicator for when the array is full
 */
template<typename T, size_t MAX_SIZE>
class CircularBuffer
{
public:
    CircularBuffer() : max_size(MAX_SIZE), head(0), tail(0), full(false) {}

    T& operator[](const size_t& idx)
    {
        if (idx < max_size)
            return buf[(head + idx) % max_size];
        else
            throw std::runtime_error("CircularBuffer: index > max_size");
    }

    friend std::ostream& operator<<(std::ostream& os, const CircularBuffer& cbuf)
    {
        os << "ptr: " << cbuf.buf << ", max_size: " << cbuf.max_size << ", size: " << cbuf.size();
        return os;
    }

    T& front() { return buf[head]; }
    T& back() { return buf[(tail - 1 + max_size) % max_size]; }
    const T& get(const size_t& idx) const { return buf[(head+idx) % max_size]; }
    const size_t size() const { return full ? max_size : (tail+max_size) % max_size; }

    void print() const
    {
        for (size_t i = 0; i < max_size; ++i)
            std::cout << get(i) << " ";
        std::cout << std::endl;
    }

    void push_back(const T& value)
    {
        if (full && head == tail)
            head = ++head % max_size;
        buf[tail] = value;
        tail = ++tail % max_size;
        if (!full && head == tail)
            full = true;
    }

    void clear()
    {
        full = false;
        head = 0;
        tail = 0;
    }

    bool isfull() {
      return full;
    }

private:
    T buf[MAX_SIZE];        // Main data buffer
    const size_t max_size;  // Size of array when full
    size_t head;            // Index to beginning of the array
    size_t tail;            // Index to the end of the array
    bool full;              // Indicator for when the array is full
};


} // namespace common

#endif // COMMON_H
