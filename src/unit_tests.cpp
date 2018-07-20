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

template<typename T>
void TEST(std::string test_name, Eigen::MatrixBase<T> mat1, Eigen::MatrixBase<T> mat2)
{
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
  Eigen::Vector3d delta;
  randomNormalMatrix(delta,dist,rng);
  Quaternion q2 = q1 + delta;
  Eigen::Vector3d delta_new = q2 - q1;
  TEST("Quaternion boxplus and boxminus operations.",delta,delta_new);


}
