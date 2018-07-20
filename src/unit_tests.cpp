#include "common_cpp/common.h"
#include <random>
#include <chrono>

int main()
{
  // initialize Gaussian random number generation
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine rng(seed);
  std::normal_distribution<double> dist(0.0,1.0);

  Eigen::Vector3d M;
  common::randomNormalMatrix(M,dist,rng);
  std::cout << M << std::endl;
}
