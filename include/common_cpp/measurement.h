#pragma once

#include <Eigen/Dense>
#include <vector>
#include <set>
#include "geometry/support.h"
#include "geometry/xform.h"


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
  ROTENC = 8
};


template<typename T>
class Imu
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int id;
  int type;
  T t;
  Eigen::Matrix<T,3,1> accel;
  Eigen::Matrix<T,3,1> gyro;

  Imu()
    : id(-1), type(IMU), t(NAN)
  {
    accel.setConstant(NAN);
    gyro.setConstant(NAN);
  }

  Imu(const int& _id, const T& _t, const Eigen::Matrix<T,3,1>& _accel, const Eigen::Matrix<T,3,1>& _gyro)
     : id(_id), type(IMU), accel(_accel), gyro(_gyro)
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
  T t;
  xform::Xform<T> transform; // transform containing position and attitude

  Mocap()
    : id(-1), type(MOCAP), t(NAN)
  {
    transform.t_.setConstant(NAN);
    transform.q_.arr_.setConstant(NAN);
  }

  Mocap(const int& _id, const T& _t, const xform::Xform<T>& _transform)
    : id(_id), type(MOCAP), t(_t)
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
  T t;
  Eigen::Matrix<T,3,1> pos; // position
  Eigen::Matrix<T,3,1> vel; // velocity

  Gps()
    : id(-1), type(GPS), t(NAN)
  {
    pos.setConstant(NAN);
    vel.setConstant(NAN);
  }

  Gps(const int& _id, const T& _t, const Eigen::Matrix<T,3,1>& _pos, const Eigen::Matrix<T,3,1>& _vel)
    : id(_id), type(GPS), t(_t), pos(_pos), vel(_vel)
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

  Feat(const int& _id, const Eigen::Matrix<T,2,1>& _pix, const Eigen::Matrix<T,3,1>& _pos)
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
  T t;
  FeatVec<T> feats; // tracked features in the image

  Image()
    : id(-1), type(IMAGE), t(NAN)
  {}

  Image(const int& _id, const T& _t, const FeatVec<T>& _feats)
    : id(_id), type(IMAGE), t(_t), feats(_feats)
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
  T t;
  T pres; // absolute pressure (Pa)

  Baro()
    : id(-1), type(BARO), t(NAN), pres(NAN)
  {}

  Baro(const int& _id, const T& _t, const T& _pres)
    : id(_id), type(BARO), t(_t), pres(_pres)
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
  T t;
  Eigen::Matrix<T,3,1> field; // magnetic field (nanotesla)

  Mag()
    : id(-1), type(MAG), t(NAN)
  {
    field.setConstant(NAN);
  }

  Mag(const int& _id, const T& _t, const Eigen::Matrix<T,3,1>& _field)
    : id(_id), type(MAG), t(_t), field(_field)
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
  T t;
  T pres; // differential pressure (Pa)

  Pitot()
    : id(-1), type(PITOT), t(NAN), pres(NAN)
  {}

  Pitot(const int& _id, const T& _t, const T& _pres)
    : id(_id), type(PITOT), t(_t), pres(_pres)
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
  T t;
  T angle; // angle (rad)

  Wvane()
    : id(-1), type(WVANE), t(NAN), angle(NAN)
  {}

  Wvane(const int& _id, const T& _t, const T& _angle)
    : id(_id), type(WVANE), t(_t), angle(_angle)
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
  T t;
  T angle; // angle (rad)

  RotEnc()
    : id(-1), type(ROTENC), t(NAN), angle(NAN)
  {}

  RotEnc(const int& _id, const T& _t, const T& _angle)
    : id(_id), type(ROTENC), t(_t), angle(_angle)
  {}
};
typedef RotEnc<float> RotEncf;
typedef RotEnc<double> RotEncd;


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

  T t;
  int type;
  Imu<T> imu;
  Mocap<T> mocap;
  Gps<T> gps;
  Image<T> image;
  Baro<T> baro;
  Mag<T> mag;
  Pitot<T> pitot;
  Wvane<T> wvane;

  Measurement(const int& _type, const T& _t, const Imu<T>& _imu)
    : type(_type), t(_t), imu(_imu)
  {}

  Measurement(const int& _type, const T& _t, const Mocap<T>& _mocap)
    : type(_type), t(_t), mocap(_mocap)
  {}

  Measurement(const int& _type, const T& _t, const Gps<T>& _gps)
    : type(_type), t(_t), gps(_gps)
  {}

  Measurement(const int& _type, const T& _t, const Image<T>& _image)
    : type(_type), t(_t), image(_image)
  {}

  Measurement(const int& _type, const T& _t, const Baro<T>& _baro)
    : type(_type), t(_t), baro(_baro)
  {}

  Measurement(const int& _type, const T& _t, const Pitot<T>& _pitot)
    : type(_type), t(_t), pitot(_pitot)
  {}

  Measurement(const int& _type, const T& _t, const Wvane<T>& _wvane)
    : type(_type), t(_t), wvane(_wvane)
  {}

  bool operator< (const Measurement<T>& other) const
  {
    return t < other.t;
  }
};
typedef Measurement<float> Measurementf;
typedef Measurement<double> Measurementd;
template<typename T2>
using Measurements = std::multiset<Measurement<T2>>;
typedef Measurements<float> Measurementsf;
typedef Measurements<double> Measurementsd;


} // namespace common
