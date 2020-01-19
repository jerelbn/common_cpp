// Quaternion unit tests.
#include <random>
#include <gtest/gtest.h>
#include "common_cpp/common.h"
#include "common_cpp/quaternion.h"

using namespace std;
using namespace Eigen;

namespace common
{

#define NUM_ITERS 10000
#define TOL 1e-6

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


static const unsigned seed = time(NULL);
static const double sqrt3 = sqrt(3.0);


TEST(Rotation, Exponential)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Matrix3d R1 = Matrix3d::Identity();
    Vector3d delta = M_PI / sqrt3 * Vector3d::Random();
    Matrix3d R2 = R1*expR(skew(delta));
    Matrix3d R3 = R2*expR(skew(-delta));
    EXPECT_MATRIX_CLOSE(R1, R3, TOL);
  }
}


TEST(Rotation, Logarithm)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Matrix3d R1 = Matrix3d::Identity();
    Vector3d delta = M_PI / sqrt3 * Vector3d::Random();
    Matrix3d R2 = R1*expR(skew(delta));
    Vector3d delta2 = vex(logR(Matrix3d(R1.transpose()*R2)));
    EXPECT_MATRIX_CLOSE(delta, delta2, TOL);
  }
}


TEST(Rotation, ExponentialTranspose)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Matrix3d R1 = Matrix3d::Identity();
    Vector3d delta = M_PI / sqrt3 * Vector3d::Random();
    Matrix3d R2 = expR(skew(delta))*R1;
    Matrix3d R3 = expR(skew(-delta))*R2;
    EXPECT_MATRIX_CLOSE(R1, R3, TOL);
  }
}


TEST(Rotation, LogarithmTranspose)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Matrix3d R1 = Matrix3d::Identity();
    Vector3d delta = M_PI / sqrt3 * Vector3d::Random();
    Matrix3d R2 = expR(skew(delta))*R1;
    Vector3d delta2 = vex(logR(Matrix3d(R1.transpose()*R2)));
    EXPECT_MATRIX_CLOSE(delta, delta2, TOL);
  }
}


} // namespace common
