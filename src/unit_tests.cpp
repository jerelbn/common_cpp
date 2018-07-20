#include "common_cpp/common.h"
#include <random>
#include <chrono>

using namespace common;

double tol = 1e-6;

std::string redText(std::string input)
{
  return "\033[31m" + input + "\033[0m";
}

std::string greenText(std::string input)
{
  return "\033[32m" + input + "\033[0m";
}

template<typename T1, typename T2>
void TEST(std::string test_name, Eigen::MatrixBase<T1>& mat1, Eigen::MatrixBase<T2>& mat2)
{
  if (mat1.norm() < 1e-6 || mat2.norm() < 1e-6)
    std::cout << redText("WARNING: Test values near zero.\n");
  if (fabs((mat1 - mat2).norm()) < tol)
  {
    std::cout << greenText("[PASSED] ") << test_name << std::endl;
  }
  else
  {
    std::cout << redText("[FAILED] ") << test_name << std::endl;
    std::cout << "\nvalue1 = \n" << mat1 << std::endl;
    std::cout << "\nvalue2 = \n" << mat2 << std::endl;
  }
}

int main()
{
  // Initialize Gaussian random number generation
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine rng(seed);
  std::normal_distribution<double> dist(0.0,1.0);

  // Check boxplus and boxminus operations
  Eigen::Vector4d vec;
  randomNormalMatrix(vec,dist,rng);
  Quaternion q1(vec);
  q1.normalize();
  Eigen::Vector3d delta;
  randomNormalMatrix(delta,dist,rng);
  Quaternion q2 = q1 + delta;
  Eigen::Vector3d delta_new = q2 - q1;
  TEST("Quaternion boxplus and boxminus operations.",delta,delta_new);

  // Check skew and vex
  delta_new = vex(skew(delta));
  TEST("Skew and vex operations.",delta,delta_new);

  // Check rotation matrix exponential
  Eigen::Matrix3d R1 = q1.R();
  Eigen::Matrix3d R2 = q2.R();
  Eigen::Vector3d deltaR = -R1.transpose()*delta;
  Eigen::Matrix3d R2_new = R1*expR(skew(deltaR));
  TEST("Rotation matrix exponential.",R2,R2_new);

  // Check rotation matrix logarithm
  Eigen::Matrix3d dR = R1*R2.transpose();
  delta_new = vex(logR(dR));
  TEST("Rotation matrix logarithm.",delta,delta_new);
}
