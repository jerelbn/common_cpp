// Quaternion unit tests.
#include <random>
#include <gtest/gtest.h>
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
static const default_random_engine rng(seed);
static const std::uniform_real_distribution<double> dist(-1.0,1.0);
static const double sqrt2 = sqrt(2.0);
static const double sqrt3 = sqrt(3.0);


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
    Quaterniond q1(u1);
    Vector2d delta = M_PI / sqrt2 * Vector2d::Random();
    Quaterniond q2 = Quaterniond::exp(q1.proj()*delta) * q1;
    Vector2d delta2 = Quaterniond::log_uvec(q2, q1);
    EXPECT_MATRIX_CLOSE(delta, delta2, TOL);
  }
}


} // namespace common
