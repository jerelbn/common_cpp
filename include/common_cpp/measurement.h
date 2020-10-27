#pragma once

#include <Eigen/Dense>
#include <vector>
#include <set>
#include "common_cpp/transform.h"


namespace common
{


enum
{
  IMU = 0,
  MOCAP = 1,
  GPS = 2,
  IMAGE = 3,
  BARO = 4,
  MAG = 5,
  PITOT = 6,
  WVANE = 7,
  ROTENC = 8,
  LRF = 9,
  LST = 10
};


template<typename T>
class Imu
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  uint32_t t_ms;
  Eigen::Matrix<T,3,1> accel;
  Eigen::Matrix<T,3,1> gyro;

  Imu()
    : id(-1), type(IMU), t_ms(0)
  {
    accel.setConstant(NAN);
    gyro.setConstant(NAN);
  }

  Imu(int _id, uint32_t _t_ms, const Eigen::Matrix<T,3,1>& _accel, const Eigen::Matrix<T,3,1>& _gyro)
     : id(_id), t_ms(_t_ms), type(IMU), accel(_accel), gyro(_gyro)
  {}

  Eigen::Matrix<T,6,1> vec() const
  {
    Eigen::Matrix<T,6,1> out;
    out << accel, gyro;
    return out;
  }
};
typedef Imu<float> Imuf;
typedef Imu<double> Imud;


template<typename T>
class Mocap
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  uint32_t t_ms;
  common::Transform<T> transform; // transform containing position and attitude

  Mocap()
    : id(-1), type(MOCAP), t_ms(0)
  {
    transform.p(Eigen::Matrix<T,3,1>::Constant(NAN));
    transform.q(Eigen::Matrix<T,4,1>::Constant(NAN));
  }

  Mocap(int _id, uint32_t _t_ms, const common::Transform<T>& _transform)
    : id(_id), type(MOCAP), t_ms(_t_ms)
  {
    transform = _transform;
  }
};
typedef Mocap<float> Mocapf;
typedef Mocap<double> Mocapd;


template<typename T>
class Gps
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  uint32_t t_ms;
  Eigen::Matrix<T,3,1> pos; // position
  Eigen::Matrix<T,3,1> vel; // velocity

  Gps()
    : id(-1), type(GPS), t_ms(0)
  {
    pos.setConstant(NAN);
    vel.setConstant(NAN);
  }

  Gps(int _id, uint32_t _t_ms, const Eigen::Matrix<T,3,1>& _pos, const Eigen::Matrix<T,3,1>& _vel)
    : id(_id), type(GPS), t_ms(_t_ms), pos(_pos), vel(_vel)
  {}

  Eigen::Matrix<T,6,1> vec() const
  {
    Eigen::Matrix<T,6,1> out;
    out << pos, vel;
    return out;
  }
};
typedef Gps<float> Gpsf;
typedef Gps<double> Gpsd;


template<typename T>
class Feat
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id; // feature id or label
  Eigen::Matrix<T,2,1> pix; // pixel position in image
  Eigen::Matrix<T,3,1> pos; // vector from camera to landmark in camera frame
  T rho; // inverse z component of pos
  T depth; // magnitude of pos

  Feat()
    : id(-1), rho(NAN), depth(NAN)
  {
    pix.setConstant(NAN);
    pos.setConstant(NAN);
  }

  Feat(int _id, const Eigen::Matrix<T,2,1>& _pix, const Eigen::Matrix<T,3,1>& _pos)
    : id(_id), pix(_pix), pos(_pos)
  {
    rho = T(1.0) / pos(2);
    depth = pos.norm();
  }
};
template<typename T2>
using FeatVec = std::vector<Feat<T2>, Eigen::aligned_allocator<Feat<T2>>>;
typedef Feat<float> Featf;
typedef Feat<double> Featd;
typedef std::vector<Featf, Eigen::aligned_allocator<Featf>> FeatVecf;
typedef std::vector<Featd, Eigen::aligned_allocator<Featd>> FeatVecd;


template<typename T>
class Image
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  uint32_t t_ms;
  FeatVec<T> feats; // tracked features in the image

  Image()
    : id(-1), type(IMAGE), t_ms(0)
  {}

  Image(int _id, uint32_t _t_ms, const FeatVec<T>& _feats)
    : id(_id), type(IMAGE), t_ms(_t_ms), feats(_feats)
  {}
};
typedef Image<float> Imagef;
typedef Image<double> Imaged;


template<typename T>
class Baro
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  uint32_t t_ms;
  T pres; // absolute pressure (Pa)

  Baro()
    : id(-1), type(BARO), t_ms(0), pres(NAN)
  {}

  Baro(int _id, uint32_t _t_ms, const T& _pres)
    : id(_id), type(BARO), t_ms(_t_ms), pres(_pres)
  {}
};
typedef Baro<float> Barof;
typedef Baro<double> Barod;


template<typename T>
class Mag
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  uint32_t t_ms;
  Eigen::Matrix<T,3,1> field; // magnetic field (nanotesla)

  Mag()
    : id(-1), type(MAG), t_ms(0)
  {
    field.setConstant(NAN);
  }

  Mag(int _id, uint32_t _t_ms, const Eigen::Matrix<T,3,1>& _field)
    : id(_id), type(MAG), t_ms(_t_ms), field(_field)
  {}
};
typedef Mag<float> Magf;
typedef Mag<double> Magd;


template<typename T>
class Pitot
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  uint32_t t_ms;
  T pres; // differential pressure (Pa)

  Pitot()
    : id(-1), type(PITOT), t_ms(0), pres(NAN)
  {}

  Pitot(int _id, uint32_t _t_ms, const T& _pres)
    : id(_id), type(PITOT), t_ms(_t_ms), pres(_pres)
  {}
};
typedef Pitot<float> Pitotf;
typedef Pitot<double> Pitotd;


template<typename T>
class Wvane
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  uint32_t t_ms;
  T angle; // angle (rad)

  Wvane()
    : id(-1), type(WVANE), t_ms(0), angle(NAN)
  {}

  Wvane(int _id, uint32_t _t_ms, const T& _angle)
    : id(_id), type(WVANE), t_ms(_t_ms), angle(_angle)
  {}
};
typedef Wvane<float> Wvanef;
typedef Wvane<double> Wvaned;


template<typename T>
class RotEnc
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  uint32_t t_ms;
  T angle; // angle (rad)

  RotEnc()
    : id(-1), type(ROTENC), t_ms(0), angle(NAN)
  {}

  RotEnc(int _id, uint32_t _t_ms, const T& _angle)
    : id(_id), type(ROTENC), t_ms(_t_ms), angle(_angle)
  {}
};
typedef RotEnc<float> RotEncf;
typedef RotEnc<double> RotEncd;


template<typename T>
class Lrf
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  uint32_t t_ms;
  T range; // range (meters)

  Lrf()
    : id(-1), type(LRF), t_ms(0)
  {
    range = NAN;
  }

  Lrf(int _id, uint32_t _t_ms, const T& _range)
    : id(_id), type(LRF), t_ms(_t_ms), range(_range)
  {}
};
typedef Lrf<float> Lrff;
typedef Lrf<double> Lrfd;


template<typename T>
class Lst
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  uint32_t t_ms;
  T az; // azimuth (radians)
  T el; // elevation (radians)
  bool valid; // invalid when spot is outside of the FOV

  Lst()
    : id(-1), type(LST), t_ms(0)
  {
    az = NAN;
    el = NAN;
    valid = false;
  }

  Lst(int _id, uint32_t _t_ms, const T& _az, const T& _el, const bool& _valid)
    : id(_id), type(LST), t_ms(_t_ms), az(_az), el(_el), valid(_valid)
  {}
};
typedef Lst<float> Lstf;
typedef Lst<double> Lstd;


template<typename T>
struct measurement_compare
{
    bool operator()(const T& lhs, const T& rhs) const
    {
        return lhs.t < rhs.t;
    }
};


template<typename T>
class Measurement
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  uint32_t t_ms;
  int type;
  Imu<T> imu;
  Mocap<T> mocap;
  Gps<T> gps;
  Image<T> image;
  Baro<T> baro;
  Mag<T> mag;
  Pitot<T> pitot;
  Wvane<T> wvane;

  Measurement(const int& _type, uint32_t _t_ms, const Imu<T>& _imu)
    : type(_type), t_ms(_t_ms), imu(_imu)
  {}

  Measurement(const int& _type, uint32_t _t_ms, const Mocap<T>& _mocap)
    : type(_type), t_ms(_t_ms), mocap(_mocap)
  {}

  Measurement(const int& _type, uint32_t _t_ms, const Gps<T>& _gps)
    : type(_type), t_ms(_t_ms), gps(_gps)
  {}

  Measurement(const int& _type, uint32_t _t_ms, const Image<T>& _image)
    : type(_type), t_ms(_t_ms), image(_image)
  {}

  Measurement(const int& _type, uint32_t _t_ms, const Baro<T>& _baro)
    : type(_type), t_ms(_t_ms), baro(_baro)
  {}

  Measurement(const int& _type, uint32_t _t_ms, const Pitot<T>& _pitot)
    : type(_type), t_ms(_t_ms), pitot(_pitot)
  {}

  Measurement(const int& _type, uint32_t _t_ms, const Wvane<T>& _wvane)
    : type(_type), t_ms(_t_ms), wvane(_wvane)
  {}

  Measurement(const int& _type, uint32_t _t_ms, const Mag<T>& _mag)
    : type(_type), t_ms(_t_ms), mag(_mag)
  {}

  bool operator< (const Measurement<T>& other) const
  {
    return t_ms < other.t_ms;
  }
};
typedef Measurement<float> Measurementf;
typedef Measurement<double> Measurementd;
template<typename T2>
using Measurements = std::multiset<Measurement<T2>>;
typedef Measurements<float> Measurementsf;
typedef Measurements<double> Measurementsd;


} // namespace common
