// Quaternion unit tests.
#include <random>
#include <gtest/gtest.h>
#include "common_cpp/quaternion.h"

using namespace std;
using namespace Eigen;


namespace common
{


static const unsigned NUM_ITERS = 100000;
static const double TOL = 1e-6;
static const unsigned seed = 0;//time(NULL);
static const double sqrt2 = sqrt(2.0);
static const double sqrt3 = sqrt(3.0);


#define EXPECT_MATRIX_CLOSE(m1, m2, tol)\
{\
  EXPECT_EQ(m1.rows(), m2.rows());\
  EXPECT_EQ(m1.cols(), m2.cols());\
  for (int i = 0; i < m1.rows(); ++i)\
  {\
    for (int j = 0; j < m1.cols(); ++j)\
    {\
      EXPECT_NEAR(m1(i,j), m2(i,j), tol);\
    }\
  }\
}


TEST(Quaternion, BoxPlus)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Quaterniond q1;
    Vector3d delta = M_PI / sqrt3 * Vector3d::Random();
    Quaterniond q2 = q1 + delta;
    Quaterniond q3 = q2 + -delta;
    EXPECT_MATRIX_CLOSE(q1.toEigen(), q3.toEigen(), TOL);
  }
}


TEST(Quaternion, BoxMinus)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Quaterniond q1;
    Vector3d delta = M_PI / sqrt3 * Vector3d::Random();
    Quaterniond q2 = q1 + delta;
    Vector3d delta2 = q2 - q1;
    EXPECT_MATRIX_CLOSE(delta, delta2, TOL);
  }
}


TEST(Quaternion, UnitVector)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Vector3d u1(0,0,1);
    Quaterniond q1 = common::Quaterniond::fromUnitVector(u1);
    Vector2d delta = M_PI / sqrt2 * Vector2d::Random();
    Quaterniond q2 = Quaterniond::exp(q1.proj()*delta) * q1;
    Vector2d delta2 = Quaterniond::logUnitVector(q2, q1);
    EXPECT_MATRIX_CLOSE(delta, delta2, TOL);
  }
}


TEST(Quaternion, Euler321)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Vector3d xyz1 = M_PI / 2 * Vector3d::Random();
    Quaterniond q1 = common::Quaterniond::fromEuler(xyz1(0), xyz1(1), xyz1(2));
    Quaterniond q2 = Quaterniond::fromYaw(xyz1(2)) * Quaterniond::fromPitch(xyz1(1)) * Quaterniond::fromRoll(xyz1(0));
    Vector3d xyz2(q2.roll(), q2.pitch(), q2.yaw());
    Vector3d xyz3(q1.roll(), q1.pitch(), q1.yaw());
    EXPECT_MATRIX_CLOSE(xyz1, xyz2, TOL);
    EXPECT_MATRIX_CLOSE(xyz1, xyz3, TOL);
  }
}


TEST(Quaternion, Euler312)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Vector3d xyz1 = M_PI / 2 * Vector3d::Random();
    Quaterniond q1 = common::Quaterniond::fromEuler(xyz1(0), xyz1(1), xyz1(2), 312);
    Quaterniond q2 = Quaterniond::fromYaw(xyz1(2)) * Quaterniond::fromRoll(xyz1(0)) * Quaterniond::fromPitch(xyz1(1));
    q2.eulerOrder(312);
    Vector3d xyz2(q2.roll(), q2.pitch(), q2.yaw());
    Vector3d xyz3(q1.roll(), q1.pitch(), q1.yaw());
    EXPECT_MATRIX_CLOSE(xyz1, xyz2, TOL);
    EXPECT_MATRIX_CLOSE(xyz1, xyz3, TOL);
  }
}


TEST(Quaternion, QuaternionFromRotationMatrix)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Quaterniond q1(Vector4d::Random().normalized());
    Quaterniond q2 = common::Quaterniond::fromRotationMatrix(q1.R());
    EXPECT_MATRIX_CLOSE(q1.toEigen(), q2.toEigen(), TOL);
  }
}


} // namespace common
