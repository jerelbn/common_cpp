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

  Quaternion(const int& order=321)
  {
    arr(0) = T(1.0);
    arr(1) = T(0.0);
    arr(2) = T(0.0);
    arr(3) = T(0.0);
    eulerOrder(order);
  }

  Quaternion(const T* ptr, const int& order=321)
  {
    arr(0) = ptr[0];
    arr(1) = ptr[1];
    arr(2) = ptr[2];
    arr(3) = ptr[3];
    if (arr(0) < T(0)) arr *= T(-1);
    eulerOrder(order);
  }

  Quaternion(const T& _w, const T& _x, const T& _y, const T& _z, const int& order=321)
  {
    arr(0) = _w;
    arr(1) = _x;
    arr(2) = _y;
    arr(3) = _z;
    if (arr(0) < T(0)) arr *= T(-1);
    eulerOrder(order);
  }

  Quaternion(const Eigen::Matrix<T,4,1>& v, const int& order=321)
  {
    arr(0) = v(0);
    arr(1) = v(1);
    arr(2) = v(2);
    arr(3) = v(3);
    if (arr(0) < T(0)) arr *= T(-1);
    eulerOrder(order);
  }

  Quaternion(const Quaternion<T>& q, const int& order=321)
  {
    arr = q.toEigen();
    if (arr(0) < T(0)) arr *= T(-1);
    eulerOrder(order);
  }

  void operator=(const Quaternion<T> &q2)
  {
    arr = q2.toEigen();
    if (arr(0) < T(0)) arr *= T(-1);
    eulerOrder(q2.eulerOrder());
  }

  Quaternion<T> operator*(const Quaternion<T> &q2) const
  {
    // if (eulerOrder() != q2.eulerOrder())
    //   std::cout << "\n\n\nWARNING: Multiplying quaterions with non-matching Euler angle rotation orders!\n\n\n";
    return Quaternion<T>(w()*q2.w() - x()*q2.x() - y()*q2.y() - z()*q2.z(),
                         w()*q2.x() + x()*q2.w() + y()*q2.z() - z()*q2.y(),
                         w()*q2.y() - x()*q2.z() + y()*q2.w() + z()*q2.x(),
                         w()*q2.z() + x()*q2.y() - y()*q2.x() + z()*q2.w(),
                         eulerOrder());
  }

  void operator*=(const Quaternion<T> &q)
  {
    *this = *this * q;
  }

  Quaternion<T> operator+(const Eigen::Matrix<T,3,1>& delta) const
  {
    Quaternion<T> q = *this * exp(delta);
    return q;
  }

  void operator+=(const Eigen::Matrix<T,3,1>& delta)
  {
    *this = *this + delta;
  }

  Eigen::Matrix<T,3,1> operator-(const Quaternion<T> &q2) const
  {
    return log(q2.inverse() * *this);
  }

  friend std::ostream& operator<<(std::ostream &os, const Quaternion<T> &q)
  {
    os << q.w() << "\n" << q.x() << "\n" << q.y() << "\n" << q.z();
    return os;
  }

  T mag() const
  {
    return sqrt(w()*w() + x()*x() + y()*y() + z()*z());
  }

  void normalize()
  {
    T m = mag();
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

  Eigen::Matrix<T,3,1> bar() const
  {
    return Eigen::Matrix<T,3,1>(x(),y(),z());
  }

  Quaternion<T> inverse() const
  {
    return Quaternion<T>(w(), -x(), -y(), -z());
  }

  Eigen::Matrix<T,3,1> rotSlow(const Eigen::Matrix<T,3,1>& v) const
  {
    Quaternion qv(T(0.0), v(0), v(1), v(2));
    Quaternion qv_new = inverse() * qv * *this;
    return Eigen::Matrix<T,3,1>(qv_new.x(), qv_new.y(), qv_new.z());
  }

  // Passive rotation (rotate coordinate frame)
  Eigen::Matrix<T,3,1> rotp(const Eigen::Matrix<T,3,1>& v) const
  {
    Eigen::Matrix<T,3,1> t = T(2.0) * v.cross(bar());
    return v + w() * t + t.cross(bar());
  }

  // Active rotation (rotate vector, transpose of passive rotation)
  Eigen::Matrix<T,3,1> rota(const Eigen::Matrix<T,3,1>& v) const
  {
    Eigen::Matrix<T,3,1> t = T(2.0) * v.cross(bar());
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

  
  /** @brief Returns x-axis rotation for 3-2-1 or 3-1-2 Euler angles.
   */
  T roll() const
  {
    if (eulerOrder() == 321)
    {
      return atan2(T(2.0) * (w()*x() + y()*z()), T(1.0) - T(2.0) * (x()*x() + y()*y()));
    }
    else if (eulerOrder() == 312)
    {
      // Hold at 90 degrees if invalid
      T val = T(2.0) * (w()*x() + y()*z());
      if (std::abs(val) > T(1.0))
        return copysign(T(1.0), val) * T(M_PI) / T(2.0);
      else
        return asin(val);
    }
  }


  /** @brief Returns y-axis rotation for 3-2-1 or 3-1-2 Euler angles.
   */
  T pitch() const
  {
    if (eulerOrder() == 321)
    {
      // Hold at 90 degrees if invalid
      const T val = T(2.0) * (w()*y() - x()*z());
      if (std::abs(val) > T(1.0))
        return copysign(T(1.0), val) * T(M_PI) / T(2.0);
      else
        return asin(val);
    }
    else if (eulerOrder() == 312)
    {
      return atan2(T(2.0) * (w()*y() - x()*z()), T(2.0) * (w()*w() + z()*z()) - T(1.0));
    }
  }


  /** @brief Returns z-axis rotation for 3-2-1 or 3-1-2 Euler angles.
   */
  T yaw() const
  {
    if (eulerOrder() == 321)
    {
      return atan2(T(2.0) * (w()*z() + x()*y()), T(1.0) - T(2.0) * (y()*y() + z()*z()));
    }
    else if (eulerOrder() == 312)
    {
      return atan2(T(2.0) * (w()*z() - x()*y()), T(2.0) * (w()*w() + y()*y()) - T(1.0));
    }
  }


  /** @brief Returns vector of Euler angles.
   */
  Eigen::Matrix<T,3,1> eulerVector() const
  {
    return Eigen::Matrix<T,3,1>(roll(), pitch(), yaw());
  }


  /** @brief Returns rotation matrix from unit quaternion.
   */
  Eigen::Matrix<T,3,3> R() const
  {
    // Pre-calculations
    T ww = w() * w();
    T wx = w() * x();
    T wy = w() * y();
    T wz = w() * z();
    T xy = x() * y();
    T xz = x() * z();
    T yz = y() * z();

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

  static Quaternion<T> fromEigen(const Eigen::Matrix<T,4,1>& q)
  {
    return Quaternion<T>(q);
  }


  /** @brief Creates quaternion from a rotation angle about the x-axis.
   @param roll rotation about x-axis
   */
  static Quaternion<T> fromRoll(const T& roll)
  {
    // Pre-calculations
    T r_2 = roll / T(2.0);

    // Output
    return Quaternion<T>(cos(r_2),
                         sin(r_2),
                         0,
                         0);
  }


  /** @brief Creates quaternion from a rotation angle about the y-axis.
   @param pitch rotation about x-axis
   */
  static Quaternion<T> fromPitch(const T& pitch)
  {
    // Pre-calculations
    T p_2 = pitch / T(2.0);

    // Output
    return Quaternion<T>(cos(p_2),
                         0,
                         sin(p_2),
                         0);
  }


  /** @brief Creates quaternion from a rotation angle about the z-axis.
   @param yaw rotation about x-axis
   */
  static Quaternion<T> fromYaw(const T& yaw)
  {
    // Pre-calculations
    T y_2 = yaw / T(2.0);

    // Output
    return Quaternion<T>(cos(y_2),
                         0,
                         0,
                         sin(y_2));
  }


  /** @brief Creates quaternion from Euler angles.
   @note Currently supports 321 and 312 rotation orders. 
   The order refers to the sequence of rotations about each axis. 
   E.g. 321 rotates the fixed coordinate frame by yaw, pitch, then roll to arrive at the body frame.
   @param roll rotation about x-axis
   @param pitch rotation about y-axis
   @param yaw rotation about z-axis
   @param order Euler angle order of rotation
   */
  static Quaternion<T> fromEuler(const T& roll, const T& pitch, const T& yaw, const int& order=321)
  {
    // Pre-calculations
    T r_2 = roll / T(2.0);
    T p_2 = pitch / T(2.0);
    T y_2 = yaw / T(2.0);
    T sr = sin(r_2);
    T sp = sin(p_2);
    T sy = sin(y_2);
    T cr = cos(r_2);
    T cp = cos(p_2);
    T cy = cos(y_2);

    // Output
    if (order == 321)
    {
      return Quaternion<T>(cr*cp*cy + sr*sp*sy,
                           sr*cp*cy - cr*sp*sy,
                           cr*sp*cy + sr*cp*sy,
                           cr*cp*sy - sr*sp*cy,
                           order);
    }
    else if (order == 312)
    {
      return Quaternion<T>(cr*cp*cy - sr*sp*sy,
                           sr*cp*cy - cr*sp*sy,
                           cr*sp*cy + sr*cp*sy,
                           cr*cp*sy + sr*sp*cy,
                           order);
    }
    else
    {
      return Quaternion<T>(order); // Throws runtime error because invalid order was supplied.
    }
  }


  /** @brief Creates quaternion from an axis of rotation and angular magnitude.
   @param axis axis of rotation
   @param angle magnitude of rotation in radians
   */
  static Quaternion<T> fromAxisAngle(const Eigen::Matrix<T,3,1>& axis, const double& angle)
  {
    return exp(angle*axis);
  }


  /** @brief Creates a unit quaternion representing rotation the z-axis to a unit vector.
   @param v unit vector input
   */
  static Quaternion<T> fromUnitVector(const Eigen::Matrix<T,3,1>& v)
  {
    Quaternion<T> q;
    Eigen::Matrix<T,3,1> ev = v.normalized(); // enforce unit length
    T angle = acos(e3.cast<T>().dot(ev));
    if (angle < T(1e-8))
    {
      q.w(T(1.0));
      q.x(T(0.0));
      q.y(T(0.0));
      q.z(T(0.0));
    }
    else
    {
      Eigen::Matrix<T,3,1> axis = (e3.cast<T>().cross(ev)).normalized();
      T half_angle = angle / T(2.0);
      Eigen::Matrix<T,3,1> qbar = axis * sin(half_angle);

      q.w(cos(half_angle));
      q.x(qbar(0));
      q.y(qbar(1));
      q.z(qbar(2));
    }

    return q;
  }

  static Quaternion<T> exp(const Eigen::Matrix<T,3,1>& delta)
  {
    T delta_norm = delta.norm();

    Quaternion<T> q;
    if (delta_norm < T(1e-8)) // avoid numerical error with approximation
    {
      q.w(T(1.0));
      q.x(delta(0) / T(2.0));
      q.y(delta(1) / T(2.0));
      q.z(delta(2) / T(2.0));
    }
    else
    {
      const T delta_norm_2 = delta_norm / T(2.0);
      const T sn = sin(delta_norm_2) / delta_norm;
      q.w(cos(delta_norm_2));
      q.x(sn * delta(0));
      q.y(sn * delta(1));
      q.z(sn * delta(2));
    }

    return q;
  }

  static Eigen::Matrix<T,3,1> log(const Quaternion<T>& q)
  {
    // get magnitude of complex portion
    Eigen::Matrix<T,3,1> qbar(q.x(), q.y(), q.z());
    T qbar_mag = qbar.norm();

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
      T delta_mag = wrapAngle(T(2.0) * atan2(qbar_mag, q.w()), T(M_PI));
      delta = delta_mag * qbar / qbar_mag;
    }

    return delta;
  }

  // q1 - q2
  static Eigen::Matrix<T,2,1> logUnitVector(const Quaternion<T>& q1, const Quaternion<T>& q2)
  {
    // get unit vectors
    Eigen::Matrix<T,3,1> e1 = q1.uvec();
    Eigen::Matrix<T,3,1> e2 = q2.uvec();

    // avoid too small of angles
    T e2T_e1 = saturate<T>(e2.dot(e1), T(1.0), T(-1.0));
    if (e2T_e1 - T(1.0) > T(-1e-16)) // same direction
      return Eigen::Matrix<T,2,1>(T(0.0), T(0.0));
    else if (e2T_e1 + T(1.0) < T(1e-16)) // opposite direction
      return Eigen::Matrix<T,2,1>(T(M_PI), T(0.0));
    else
    {
      // compute axis angle difference
      Eigen::Matrix<T,3,1> e2_x_e1 = e2.cross(e1);
      Eigen::Matrix<T,3,1> s = acos(e2T_e1) * e2_x_e1.normalized();

      // place error on first vector's tangent space
      return q1.proj().transpose() * s;
    }
  }

  static Quaternion<T> boxPlusUnitVector(const Quaternion<T>& q, const Eigen::Matrix<T,2,1>& delta)
  {
    return Quaternion<T>::exp(q.proj() * delta) * q;
  }

  // derivative of quaternion exponential map
  static Eigen::Matrix<T,3,3> expDerivative(const Eigen::Matrix<T,3,1>& delta)
  {
    T dmag = delta.norm();
    Eigen::Matrix<T,3,3> delta_x = skew(delta);
    if (dmag < T(1e-8))
      return I_3x3 - T(0.5) * delta_x;
    else
    {
      T dmag2 = dmag * dmag;
      return I_3x3 - (T(1.0) - cos(dmag)) / dmag2 * delta_x +
             (dmag - sin(dmag)) / (dmag2 * dmag) * delta_x * delta_x;
    }
  }

  static Quaternion<T> fromRotationMatrix(const Eigen::Matrix<T,3,3>& R)
  {
    // To improve numerical accuracy, choose solution with largest value in the root
    T root_w = T(1.0) + R(0,0) + R(1,1) + R(2,2);
    T root_x = T(1.0) + R(0,0) - R(1,1) - R(2,2);
    T root_y = T(1.0) - R(0,0) + R(1,1) - R(2,2);
    T root_z = T(1.0) - R(0,0) - R(1,1) + R(2,2);

    T w, x, y, z;
    if (root_w >= root_x &&
        root_w >= root_y &&
        root_w >= root_z)
    {
      w = T(0.5) * sqrt(root_w);
      x = T(0.25) * (R(1,2) - R(2,1)) / w;
      y = T(0.25) * (R(2,0) - R(0,2)) / w;
      z = T(0.25) * (R(0,1) - R(1,0)) / w;
    }
    if (root_x >= root_w &&
        root_x >= root_y &&
        root_x >= root_z)
    {
      x = T(0.5) * sqrt(root_x);
      w = T(0.25) * (R(1,2) - R(2,1)) / x;
      y = T(0.25) * (R(0,1) + R(1,0)) / x;
      z = T(0.25) * (R(0,2) + R(2,0)) / x;
    }
    if (root_y >= root_w &&
        root_y >= root_x &&
        root_y >= root_z)
    {
      y = T(0.5) * sqrt(root_y);
      w = T(0.25) * (R(2,0) - R(0,2)) / y;
      x = T(0.25) * (R(1,0) + R(0,1)) / y;
      z = T(0.25) * (R(1,2) + R(2,1)) / y;
    }
    if (root_z >= root_w &&
        root_z >= root_x &&
        root_z >= root_y)
    {
      z = T(0.5) * sqrt(root_z);
      w = T(0.25) * (R(0,1) - R(1,0)) / z;
      x = T(0.25) * (R(2,0) + R(0,2)) / z;
      y = T(0.25) * (R(2,1) + R(1,2)) / z;
    }

    return Quaternion<T>(w,x,y,z);
  }


  const T& w() const { return arr(0); }
  const T& x() const { return arr(1); }
  const T& y() const { return arr(2); }
  const T& z() const { return arr(3); }
  const int& eulerOrder() const { return euler_order; }
  const Eigen::Matrix<T,4,1> toEigen() const { return arr; }
  void w(const T& _w) { arr(0) = _w; }
  void x(const T& _x) { arr(1) = _x; }
  void y(const T& _y) { arr(2) = _y; }
  void z(const T& _z) { arr(3) = _z; }
  void eulerOrder(const int& order)
  {
    if (order == 321 || order == 312)
    {
      euler_order = order;
    }
    else
    {
      std::stringstream ss;
      ss << "Invalid Euler angle order in " << __FILE__ \
         << " line " << __LINE__ << ": " << "Choose 321 or 312." << std::endl;
      throw std::runtime_error(ss.str());
    }
  }
  T* data() { return arr.data(); }

private:

  Eigen::Matrix<T,4,1> arr;
  int euler_order = 321;

};

typedef Quaternion<float> Quaternionf;
typedef Quaternion<double> Quaterniond;

} // namespace common