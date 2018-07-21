// Name: Jerel Nielsen
// Date: 15 June 2017
// Desc: Container for common functions.

#include "common_cpp/common.h"

namespace common
{


Quaternion::Quaternion() 
{
  w = 1;
  x = 0;
  y = 0;
  z = 0;
}

Quaternion::~Quaternion() {}

Quaternion::Quaternion(double _w, double _x, double _y, double _z)
{
  w = _w;
  x = _x;
  y = _y;
  z = _z;
}

Quaternion::Quaternion(double roll, double pitch, double yaw)
{
  w = cos(roll/2)*cos(pitch/2)*cos(yaw/2) + sin(roll/2)*sin(pitch/2)*sin(yaw/2);
  x = sin(roll/2)*cos(pitch/2)*cos(yaw/2) - cos(roll/2)*sin(pitch/2)*sin(yaw/2);
  y = cos(roll/2)*sin(pitch/2)*cos(yaw/2) + sin(roll/2)*cos(pitch/2)*sin(yaw/2);
  z = cos(roll/2)*cos(pitch/2)*sin(yaw/2) - sin(roll/2)*sin(pitch/2)*cos(yaw/2);
}

Quaternion::Quaternion(const Eigen::Vector4d &v)
{
  w = v(0);
  x = v(1);
  y = v(2);
  z = v(3);
}

Quaternion::Quaternion(Eigen::Vector3d &fz)
{
  // convert to axis-angle representation
  Eigen::Vector3d ez(0, 0, 1);
  fz = fz/fz.norm(); // make sure unit vector
  double theta = acos(fz.dot(ez));

  if (theta < 1e-8)
  {
    w = 1;
    x = 0;
    y = 0;
    z = 0;
  }
  else
  {
    Eigen::Vector3d fzcross_ez = fz.cross(ez);
    Eigen::Vector3d iaa = fzcross_ez/(fzcross_ez.norm());

    // get complex portion of quaternion
    Eigen::Vector3d qv = iaa * sin(theta/2);

    w = cos(theta/2);
    x = qv(0);
    y = qv(1);
    z = qv(2);
  }
}

// overload multiply operator for simple quaternion multiplication
Quaternion Quaternion::operator*(const Quaternion &q2)
{
  double qw = w*q2.w - x*q2.x - y*q2.y - z*q2.z;
  double qx = w*q2.x + x*q2.w + y*q2.z - z*q2.y;
  double qy = w*q2.y - x*q2.z + y*q2.w + z*q2.x;
  double qz = w*q2.z + x*q2.y - y*q2.x + z*q2.w;
  
  return Quaternion(qw, qx, qy, qz);
}

// overload addition operator as boxplus for a quaternion and a 3-vector
Quaternion Quaternion::operator+(const Eigen::Vector3d &delta)
{
  return (*this)*exp(delta);
}

// overload minus operator as boxminus for two quaternions
Eigen::Vector3d Quaternion::operator-(Quaternion &q2)
{
  return log(q2.inv()*(*this));
}

// overload stream operator for simple quaternion displaying
std::ostream& operator<<(std::ostream &os, const Quaternion &q)
{
  os << q.w << ", " << q.x << ", " << q.y << ", " << q.z;
  return os;
}

// quaternion inverse
Quaternion Quaternion::inv()
{
  Quaternion q;
  q.w = w;
  q.x = -x;
  q.y = -y;
  q.z = -z;

  return q;
}

// quaternion norm
double Quaternion::mag()
{
  return sqrt(w*w + x*x + y*y + z*z);
}

// quaternion normalization
void Quaternion::normalize()
{
  double mag = this->mag();
  w /= mag;
  x /= mag;
  y /= mag;
  z /= mag;

  // ensure positive scalar portion
  if (w < 0)
    w *= -1;
}

// conversion from quaternion to roll angle
double Quaternion::roll()
{
  return atan2(2*(w*x + y*z), 1 - 2*(x*x + y*y));
}

// conversion from quaternion to pitch angle
double Quaternion::pitch()
{
  double val = 2*(w*y - x*z);

  // hold at 90 degrees if invalid
  if (fabs(val) > 1)
    return copysign(1, val)*M_PI/2;
  else
    return asin(val);
}

// conversion from quaternion to yaw angle
double Quaternion::yaw()
{
  return atan2(2*(w*z + x*y), 1 - 2*(y*y + z*z));
}

// get vector of Euler angles
Eigen::Vector3d Quaternion::euler()
{
  return Eigen::Vector3d(this->roll(),this->pitch(),this->yaw());
}

// get complex portion of quaternion in Eigen
Eigen::Vector3d Quaternion::bar()
{
  Eigen::Vector3d v;
  v << x, y, z;
  return v;
}

// convert Quaternion to Eigen vector
Eigen::Vector4d Quaternion::toEigen()
{
  Eigen::Vector4d v;
  v << w, x, y, z;
  return v;
}

// convert Eigen Vector to Quaternion
void Quaternion::fromEigen(const Eigen::Vector4d &q)
{
  w = q(0);
  x = q(1);
  y = q(2);
  z = q(3);
}

// create rotation matrix from quaternion
Eigen::Matrix3d Quaternion::R()
{
  Eigen::Matrix3d R;
  R <<  2*w*w + 2*x*x - 1,      2*w*z + 2*x*y,     -2*w*y + 2*x*z,
           -2*w*z + 2*x*y,  2*w*w + 2*y*y - 1,      2*w*x + 2*y*z,
            2*w*y + 2*x*z,     -2*w*x + 2*y*z,  2*w*w + 2*z*z - 1;
  return R;
}

// rotate a 3-vector directly the original way
Eigen::Vector3d Quaternion::rotSlow(Eigen::Vector3d v)
{
  Quaternion qv = Quaternion(0, v(0), v(1), v(2));
  Quaternion qv_new = this->inv() * qv * *this;
  return Eigen::Vector3d(qv_new.x, qv_new.y, qv_new.z);
}

// rotate a 3-vector directly the fast way
Eigen::Vector3d Quaternion::rot(Eigen::Vector3d v)
{
  Eigen::Vector3d t = 2 * skew(v) * this->bar();
  return v + this->w * t + skew(t) * this->bar();
}

// compute the unit vector in the camera frame given its quaternion
Eigen::Vector3d Quaternion::uvec()
{
  Eigen::Vector3d ez;
  ez << 0, 0, 1;
  return R() * ez;
}

// projection matrix of unit vector onto its tangent space
Eigen::MatrixXd Quaternion::proj()
{
  Eigen::MatrixXd E(3,2);
  E << 1, 0,
       0, 1,
       0, 0;
  return R() * E;
}

// exponential map to unit quaternion
Quaternion Quaternion::exp(const Eigen::Vector3d &delta)
{
  double delta_norm = delta.norm();

  Quaternion q;
  if (delta_norm < 1e-8) // avoid numerical error with approximation
  {
    q.w = 1;
    q.x = delta(0)/2;
    q.y = delta(1)/2;
    q.z = delta(2)/2;
  }
  else
  {
    double sn = sin(delta_norm/2)/delta_norm;
    q.w = cos(delta_norm/2);
    q.x = sn*delta(0);
    q.y = sn*delta(1);
    q.z = sn*delta(2);
  }

  return q;
}

// unit quaternion logarithmic map to vector
Eigen::Vector3d Quaternion::log(const Quaternion &q)
{
  // get magnitude of complex portion
  Eigen::Vector3d qbar(q.x, q.y, q.z);
  double qbar_mag = qbar.norm();

  // avoid numerical error with approximation
  Eigen::Vector3d delta;
  if (qbar_mag < 1e-6)
    q.w >= 0 ? delta = qbar : delta = -qbar;
  else
    delta = 2. * atan2(qbar_mag,q.w) * qbar / qbar_mag;

  return delta;
}


} // namespace common


