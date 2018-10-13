#include "common_cpp/common.h"
#include <random>
#include <chrono>

using namespace std;
using namespace common;



int main()
{
  // Error tolerance and number of tests
  double tol = 1e-6;
  int N = 1000;
  srand((unsigned)time(NULL));

  // Initialize Gaussian random number generation
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine rng(seed);
  srand(seed);
  std::normal_distribution<double> dist(0.0,0.1);

  for (int i = 0; i < N; ++i)
  {
    Eigen::Vector3d p1; p1.setRandom();
    Quaterniond q1(dist,rng);
    Eigen::Vector3d delta_p, delta_q;
    delta_p.setRandom();
    delta_q.setRandom();
    delta_p *= 0.5;
    delta_q *= 0.5;
    if (delta_q.norm() > M_PI) continue;

    // Check boxplus and boxminus operations
    Quaterniond q2 = q1 + delta_q;
    Eigen::Vector3d delta_q_new = q2 - q1;
    if (!TEST("Quaternion boxplus and boxminus operations.",tol,delta_q,delta_q_new)) break;

    // Check rotation matrix exponential
    Eigen::Matrix3d R1 = q1.R();
    Eigen::Matrix3d R2 = q2.R();
    Eigen::Vector3d deltaR = -R1.transpose()*delta_q;
    Eigen::Matrix3d R2_new = R1*expR(skew(deltaR));
    if (!TEST("Rotation matrix exponential.",tol,R2,R2_new)) break;

    // Check rotation matrix logarithm
    Eigen::Matrix3d dR = R1*R2.transpose();
    delta_q_new = vex(logR(dR));
    if (!TEST("Rotation matrix logarithm.",tol,delta_q,delta_q_new)) break;

    // Check rotation matrix transpose exponential
    R1 = q1.R().transpose();
    R2 = q2.R().transpose();
    R2_new = R1*expR(skew(delta_q));
    if (!TEST("Rotation matrix transpose exponential.",tol,R2,R2_new)) break;

    // Check rotation matrix transpose logarithm
    dR = R1.transpose()*R2;
    delta_q_new = vex(logR(dR));
    if (!TEST("Rotation matrix transpose logarithm.",tol,delta_q,delta_q_new)) break;

    // Check transform boxplus and boxminus operators
    Transformd t1(p1, q1);
    Eigen::Matrix<double,6,1> delta_t;
    delta_t.segment<3>(0) = delta_p;
    delta_t.segment<3>(3) = delta_q;
    Transformd t2 = t1 + delta_t;
    Eigen::Matrix<double,6,1> delta_t2 = t2 - t1;
    if (!TEST("Transform boxplus and boxminus operators.",tol,delta_t,delta_t2)) break;
  }
}
