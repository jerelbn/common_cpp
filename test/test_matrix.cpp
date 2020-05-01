// Quaternion unit tests.
#include <random>
#include <gtest/gtest.h>
#include "common_cpp/matrix.h"

using namespace std;


namespace common
{


static const unsigned NUM_ITERS = 10;
static const double TOL = 1e-6;
static const unsigned seed = time(NULL);
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


TEST(Matrix, Equals)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    Matrix3d m1 = Matrix3d::random();
    Matrix3d m2 = Matrix3d::random();
    m2 = m1;
    EXPECT_MATRIX_CLOSE(m1, m2, TOL);
  }
}


TEST(Matrix, Multiplication)
{
  srand(seed);
  for (size_t iter = 0; iter < NUM_ITERS; ++iter)
  {
    double m1_[] = {1,2,3,4,5,6,7,8,9,10,11,12};
    double m2_[] = {13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
    double m4_[] = {230,240,250,260,270,558,584,610,636,662,886,928,970,1012,1054};
    Matrix<double,3,4> m1(m1_);
    Matrix<double,4,5> m2(m2_);
    Matrix<double,3,5> m3 = m1 * m2;
    Matrix<double,3,5> m4(m4_);
    EXPECT_MATRIX_CLOSE(m3, m4, TOL);
  }
}


} // namespace common
