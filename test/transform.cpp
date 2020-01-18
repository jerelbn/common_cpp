// Quaternion unit tests.
#include <random>
#include <gtest/gtest.h>
#include "common_cpp/transform.h"

using namespace std;
using namespace Eigen;

typedef Matrix<double, 6, 1> Vector6d;


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
static const double sqrt3 = sqrt(3.0);


TEST(Transform, BoxPlus)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Transformd x1;
    Vector6d delta;
    delta.head<3>() = 10.0 * Vector3d::Random();
    delta.tail<3>() = M_PI / sqrt3 * Vector3d::Random();
    Transformd x2 = x1 + delta;
    Transformd x3 = x2 + -delta;
    EXPECT_MATRIX_CLOSE(x1.toEigen(), x3.toEigen(), TOL);
  }
}


TEST(Transform, BoxMinus)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Transformd x1;
    Vector6d delta;
    delta.head<3>() = 10.0 * Vector3d::Random();
    delta.tail<3>() = M_PI / sqrt3 * Vector3d::Random();
    Transformd x2 = x1 + delta;
    Vector6d delta2 = x2 - x1;
    EXPECT_MATRIX_CLOSE(delta, delta2, TOL);
  }
}


} // namespace common
