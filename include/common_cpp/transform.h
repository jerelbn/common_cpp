#pragma once

#include <eigen3/Eigen/Eigen>
#include "common.h"
#include "quaternion.h"


namespace common
{

// This class operates on a translation p^a_b/a and rotation R^b_a
template<typename T = double>
class Transform
{

enum
{
  PX, PY, PZ, QW, QX, QY, QZ, T_SIZE
};

enum
{
  TX, TY, TZ, WX, WY, WZ, DT_SIZE
};

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Transform()
  {
    arr.setZero();
    arr(QW) = T(1.0);
  }

  Transform(const T* ptr)
  {
    arr(0) = ptr[0];
    arr(1) = ptr[1];
    arr(2) = ptr[2];
    arr(3) = ptr[3];
    arr(4) = ptr[4];
    arr(5) = ptr[5];
    arr(6) = ptr[6];
  }

  Transform(const Eigen::Matrix<T,3,1>& p, const Quaternion<T>& q)
  {
    setP(p);
    setQ(q);
  }

  Transform(const Eigen::Matrix<T,3,1>& p, const Eigen::Matrix<T,4,1>& q)
  {
    setP(p);
    setQ(q.normalized());
  }

  Transform(const T& x, const T& y, const T& z, const T& roll, const T& pitch, const T& yaw)
  {
    setPX(x);
    setPY(y);
    setPZ(z);
    setQ(Quaternion<T>(roll, pitch, yaw));
  }

  Transform(const Transform<T>& T_)
  {
    setP(T_.p());
    setQ(T_.q());
  }

  Transform(const Eigen::Matrix<T,T_SIZE,1>& t)
  {
    setT(t);
  }

  Transform<T> operator*(const Transform<T>& t2) const
  {
    Transform<T> t;
    t.setP(p() + q().inv().rot(t2.p()));
    t.setQ(q() * t2.q());
    return t;
  }

  Transform<T> operator+(const Eigen::Matrix<T,DT_SIZE,1>& delta) const
  {
    return *this * exp(delta);
  }

  void operator+=(const Eigen::Matrix<T,6,1>& delta)
  {
    *this = *this + delta;
  }

  Eigen::Matrix<T,DT_SIZE,1> operator-(const Transform<T>& t1) const
  {
    return log(t1.inv() * *this);
  }

  template<typename T2>
  Transform<T2> cast() const
  {
    return Transform<T2>(arr.template cast<T2>());
  }

  Transform<T> inv() const
  {
    return Transform<T>(-q().rot(p()), q().inv());
  }

  Eigen::Matrix<T,3,1> transform(const Eigen::Matrix<T,3,1>& v) const
  {
    return q().rot(v - p());
  }

  static Transform<T> exp(const Eigen::Matrix<T,DT_SIZE,1>& delta)
  {
    Eigen::Matrix<T,3,1> delta_t = delta.template segment<3>(TX);
    Eigen::Matrix<T,3,1> delta_q = delta.template segment<3>(WX);
    T theta = delta_q.norm();

    Transform<T> t;
    if (theta < T(1e-6)) // avoid numerical error with approximation
    {
      Eigen::Matrix<T,3,1> delta_q_2 = delta_q / T(2.0);
      t.setP(delta_t);
      t.setQ(Eigen::Matrix<T,4,1>(T(1.0), delta_q_2(0), delta_q_2(1), delta_q_2(2)));
    }
    else
    {
      T theta2 = theta * theta;
      T theta3 = theta * theta2;
      Eigen::Matrix<T,3,3> delta_q_skew = skew(delta_q);
      Eigen::Matrix<T,3,3> delta_q_skew2 = delta_q_skew * delta_q_skew;
      Eigen::Matrix<T,3,3> V = I_3x3.cast<T>() + (T(1.0) - cos(theta))/theta2 *
                               delta_q_skew + (theta - sin(theta))/theta3 * delta_q_skew2;
      t.setP(V * delta_t);
      t.setQ(Quaternion<T>::exp(delta_q));
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

  const Eigen::Matrix<T,3,1> p() const { return arr.template segment<3>(PX); }
  const Quaternion<T> q() const { return Quaternion<T>(arr(QW), arr(QX), arr(QY), arr(QZ)); }
  const Eigen::Matrix<T,T_SIZE,1> toEigen() const { return arr; }
  T* data() { return arr.data(); }

private:

  Eigen::Matrix<T,T_SIZE,1> arr;

};

typedef Transform<float> Transformf;
typedef Transform<double> Transformd;

} // namespace common