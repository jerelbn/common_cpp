#include "common_cpp/common.h"
#include <random>
#include <chrono>

using namespace common;



int main()
{
  // Error tolerance and number of tests
  double tol = 1e-6;
  int N = 1000;

  // Initialize Gaussian random number generation
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine rng(seed);
  std::normal_distribution<double> dist(0.0,1.0);

  for (int i = 0; i < N; ++i)
  {
    Eigen::Vector4d vec;
    randomNormalMatrix(vec,dist,rng);
    Quaternion q1(vec);
    q1.normalize();
    Eigen::Vector3d delta;
    randomNormalMatrix(delta,dist,rng);
    if (delta.norm() > M_PI) continue;

    // Check boxplus and boxminus operations
    Quaternion q2 = q1 + delta;
    Eigen::Vector3d delta_new = q2 - q1;
    if (!TEST("Quaternion boxplus and boxminus operations.",tol,delta,delta_new)) break;

    // Check rotation matrix exponential
    Eigen::Matrix3d R1 = q1.R();
    Eigen::Matrix3d R2 = q2.R();
    Eigen::Vector3d deltaR = -R1.transpose()*delta;
    Eigen::Matrix3d R2_new = R1*expR(skew(deltaR));
    if (!TEST("Rotation matrix exponential.",tol,R2,R2_new)) break;

    // Check rotation matrix logarithm
    Eigen::Matrix3d dR = R1*R2.transpose();
    delta_new = vex(logR(dR));
    if (!TEST("Rotation matrix logarithm.",tol,delta,delta_new)) break;

    // Check rotation matrix transpose exponential
    R1 = q1.R().transpose();
    R2 = q2.R().transpose();
    R2_new = R1*expR(skew(delta));
    if (!TEST("Rotation matrix transpose exponential.",tol,R2,R2_new)) break;

    // Check rotation matrix transpose logarithm
    dR = R1.transpose()*R2;
    delta_new = vex(logR(dR));
    if (!TEST("Rotation matrix transpose logarithm.",tol,delta,delta_new)) break;
  }
}
