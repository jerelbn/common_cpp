#include "common_cpp/common.h"
#include <random>
#include <chrono>

using namespace common;



int main()
{
  // Error tolerance 
  double tol = 1e-6;

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
  TEST("Quaternion boxplus and boxminus operations.",tol,delta,delta_new);

  // Check skew and vex
  delta_new = vex(skew(delta));
  TEST("Skew and vex operations.",tol,delta,delta_new);

  // Check rotation matrix exponential
  Eigen::Matrix3d R1 = q1.R();
  Eigen::Matrix3d R2 = q2.R();
  Eigen::Vector3d deltaR = -R1.transpose()*delta;
  Eigen::Matrix3d R2_new = R1*expR(skew(deltaR));
  TEST("Rotation matrix exponential.",tol,R2,R2_new);

  // Check rotation matrix logarithm
  Eigen::Matrix3d dR = R1*R2.transpose();
  delta_new = vex(logR(dR));
  TEST("Rotation matrix logarithm.",tol,delta,delta_new);
}
