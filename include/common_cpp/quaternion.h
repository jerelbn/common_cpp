#pragma once

#include <eigen3/Eigen/Eigen>
#include "common.h"



namespace common
{

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

  Quaternion(const T* ptr)
  {
    arr(0) = ptr[0] < T(0) ? -ptr[0] : ptr[0];
    arr(1) = ptr[1];
    arr(2) = ptr[2];
    arr(3) = ptr[3];
  }

  Quaternion(const T& _w, const T& _x, const T& _y, const T& _z)
  {
    arr(0) = _w < T(0) ? -_w : _w;
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
    arr(0) = v(0) < T(0) ? -v(0) : v(0);
    arr(1) = v(1);
    arr(2) = v(2);
    arr(3) = v(3);
  }

  // create quaternion from a unit vector
  Quaternion(Eigen::Matrix<T,3,1>& fz)
  {
    // convert to axis-angle representation
    fz.normalize(); // enforce unit length
    const T theta = acos(e3.cast<T>().dot(fz));
    if (theta < T(1e-6))
    {
      arr(0) = T(1.0);
      arr(1) = T(0.0);
      arr(2) = T(0.0);
      arr(3) = T(0.0);
    }
    else
    {
      const Eigen::Matrix<T,3,1> iaa = (e3.cast<T>().cross(fz)).normalized();
      const T theta_2 = theta / T(2.0);
      const Eigen::Matrix<T,3,1> qv = iaa * sin(theta_2);

      arr(0) = cos(theta_2);
      arr(1) = qv(0);
      arr(2) = qv(1);
      arr(3) = qv(2);
    }
  }

  // initialize random unit quaternion
  Quaternion(std::uniform_real_distribution<T>& dist, std::default_random_engine& rng)
  {
    Quaternion<T> q(dist(rng), dist(rng), dist(rng), dist(rng));
    q.normalize();
    if (q.w() < T(0)) q.setW(-q.w());
    arr(0) = q.w();
    arr(1) = q.x();
    arr(2) = q.y();
    arr(3) = q.z();
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

  Quaternion<T> operator+(const Eigen::Matrix<T,3,1>& delta) const
  {
    Quaternion<T> q = *this * exp(delta);
    if (q.w() < T(0)) q.setW(-q.w());
    return q;
  }

  void operator+=(const Eigen::Matrix<T,3,1>& delta)
  {
    *this = *this + delta;
  }

  // overload minus operator as boxminus for two quaternions
  Eigen::Matrix<T,3,1> operator-(const Quaternion<T> &q2) const
  {
    return log(q2.inverse() * *this);
  }

  friend std::ostream& operator<<(std::ostream &os, const Quaternion<T> &q)
  {
    os << q.w() << "\n" << q.x() << "\n" << q.y() << "\n" << q.z();
    return os;
  }

  static Quaternion<T> from_euler(const T& roll, const T& pitch, const T& yaw)
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
    return Quaternion<T>(cr*cp*cy + sr*sp*sy,
                         sr*cp*cy - cr*sp*sy,
                         cr*sp*cy + sr*cp*sy,
                         cr*cp*sy - sr*sp*cy);
  }

  static Quaternion<T> from_axis_angle(const Eigen::Matrix<T,3,1>& axis, const double& angle)
  {
    return exp(angle*axis);
  }

  void normalize()
  {
    const T m = mag();
    arr(0) /= m;
    arr(1) /= m;
    arr(2) /= m;
    arr(3) /= m;
  }

  Quaternion<T> normalized() const
  {
    return Quaternion<T>(Eigen::Matrix<T,4,1>(arr.normalized()));
  }

  template<typename T2>
  Quaternion<T2> cast() const
  {
    return Quaternion<T2>(arr.template cast<T2>());
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

  Quaternion<T> inverse() const
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
    const Quaternion qv_new = inverse() * qv * *this;
    return Eigen::Matrix<T,3,1>(qv_new.x(), qv_new.y(), qv_new.z());
  }

  // Passive rotation (rotate coordinate frame)
  Eigen::Matrix<T,3,1> rotp(const Eigen::Matrix<T,3,1>& v) const
  {
    const Eigen::Matrix<T,3,1> t = T(2.0) * v.cross(bar());
    return v + w() * t + t.cross(bar());
  }

  // Active rotation (rotate vector, transpose of passive rotation)
  Eigen::Matrix<T,3,1> rota(const Eigen::Matrix<T,3,1>& v) const
  {
    const Eigen::Matrix<T,3,1> t = T(2.0) * v.cross(bar());
    return v - w() * t + t.cross(bar());
  }

  Eigen::Matrix<T,3,1> uvec() const
  {
    return rota(e3.cast<T>());
  }

  Eigen::Matrix<T,3,2> proj() const
  {
    return inverse().R() * I_2x3.cast<T>().transpose();
  }

  static Quaternion<T> exp(const Eigen::Matrix<T,3,1>& delta)
  {
    const T delta_norm = delta.norm();

    Quaternion<T> q;
    if (delta_norm < T(1e-8)) // avoid numerical error with approximation
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
    if (qbar_mag < T(1e-8))
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
    const T e2T_e1 = saturate<T>(e2.dot(e1), T(1.0), T(-1.0));
    if (e2T_e1 - T(1.0) > T(-1e-16)) // same direction
      return Eigen::Matrix<T,2,1>(T(0.0), T(0.0));
    else if (e2T_e1 + T(1.0) < T(1e-16)) // opposite direction
      return Eigen::Matrix<T,2,1>(T(M_PI), T(0.0));
    else
    {
      // compute axis angle difference
      const Eigen::Matrix<T,3,1> e2_x_e1 = e2.cross(e1);
      const Eigen::Matrix<T,3,1> s = acos(e2T_e1) * e2_x_e1.normalized();

      // place error on first vector's tangent space
      return q1.proj().transpose() * s;
    }
  }

  static Quaternion<T> boxplus_uvec(const Quaternion<T>& q, const Eigen::Matrix<T,2,1>& delta)
  {
    return Quaternion<T>::exp(q.proj() * delta) * q;
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

  const T& w() const { return arr(0); }
  const T& x() const { return arr(1); }
  const T& y() const { return arr(2); }
  const T& z() const { return arr(3); }
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

} // namespace common