// Name: Jerel Nielsen
// Date: 15 June 2017
// Desc: Container for common functions.

#ifndef COMMON_H
#define COMMON_H


#include <iostream>
#include <yaml-cpp/yaml.h>
#include <eigen3/Eigen/Eigen>

namespace common
{


template<typename T>
static const T gravity = 9.80665;

class Quaternion
{

public:

  Quaternion();
  ~Quaternion();
  Quaternion(double _w, double _x, double _y, double _z);
  Quaternion(double roll, double pitch, double yaw);
  Quaternion(Eigen::Vector4d v);
  Quaternion(Eigen::Vector3d fz);

  double w;
  double x;
  double y;
  double z;

private:

  Quaternion operator*(const Quaternion &q2);
  friend std::ostream& operator<<(std::ostream &os, const Quaternion &q);
  Quaternion inv();
  double mag();
  void normalize();
  double roll();
  double pitch();
  double yaw();
  void convertFromEigen(const Eigen::Vector4d q);
  Eigen::Vector4d convertToEigen();
  Eigen::Vector3d bar();
  Eigen::Matrix3d rot();
  Eigen::Vector3d rotateVectorSlow(Eigen::Vector3d v);
  Eigen::Vector3d rotateVector(Eigen::Vector3d v);
  Eigen::Vector3d unitVector();
  Eigen::MatrixXd projection();
  Quaternion exp(const Eigen::Vector3d delta);
  Eigen::Vector3d log(const Quaternion q);

};

Eigen::VectorXd rk5(Eigen::VectorXd state, Eigen::VectorXd input, std::function<Eigen::VectorXd(Eigen::VectorXd, Eigen::VectorXd)> ode, double h);
Eigen::Vector3d log_R(const Eigen::Matrix3d R);
common::Quaternion vec2quat(const Eigen::Vector3d v);
Eigen::Vector3d vex(const Eigen::Matrix3d mat);
Eigen::Matrix3d skew(const Eigen::Vector3d vec);
Eigen::Matrix3d R_v2_to_b(double phi);
Eigen::Matrix3d R_v1_to_v2(double theta);
Eigen::Matrix3d R_v_to_v1(double psi);
Eigen::Matrix3d R_v_to_b(double phi, double theta, double psi);
Eigen::Matrix3d R_cb2c();

// Removes elements from Eigen vectors.
// idx: index of first element to remove
// num: number of elements to remove including element located at idx
template<typename T>
void removeElements(Eigen::Matrix<T,-1,1,0,-1,1>& vector, unsigned idx, unsigned num)
{
  unsigned new_elements = vector.size()-num; // number of elements in resulting vector
  if (idx < new_elements)
    vector.segment(idx,new_elements-idx) = vector.segment(idx+num,new_elements-idx);
  else
    std::cout << "ERROR: cannot remove requested elements in function 'removeElements()'!\n";
  vector.conservativeResize(new_elements);
}

// Removes rows from Eigen matrices.
// idx: index of first row to remove
// num: number of rows to remove including row located at idx
template<typename T>
void removeRows(Eigen::Matrix<T,-1,-1,0,-1,-1>& matrix, unsigned idx, unsigned num)
{
  unsigned new_rows = matrix.rows()-num; // number of rows in resulting matrix
  if (idx < new_rows)
    matrix.block(idx,0,new_rows-idx,matrix.cols()) = matrix.block(idx+num,0,new_rows-idx,matrix.cols());
  else
    std::cout << "ERROR: cannot remove requested rows in function 'removeRows()'!\n";
  matrix.conservativeResize(new_rows,matrix.cols());
}

// Removes columns from Eigen matrices.
// idx: index of first column to remove
// num: number of columns to remove including column located at idx
template<typename T>
void removeCols(Eigen::Matrix<T,-1,-1,0,-1,-1>& matrix, unsigned idx, unsigned num)
{
  unsigned new_cols = matrix.cols()-num; // number of columns in resulting matrix
  if (idx < new_cols)
    matrix.block(0,idx,matrix.rows(),new_cols-idx) = matrix.block(0,idx+num,matrix.rows(),new_cols-idx);
  else
    std::cout << "ERROR: cannot remove requested columns in function 'removeCols()'!\n";
  matrix.conservativeResize(matrix.rows(),new_cols);
}

// Loads scalar parameters from a .yaml file
// Author: James Jackson
template <typename T>
bool get_yaml_node(const std::string key, const std::string filename, T& val, bool print_error = true) 
{
  YAML::Node node = YAML::LoadFile(filename);
  if (node[key])
  {
    val = node[key].as<T>();
    return true;
  }
  else
  {
    if (print_error)
    {
      printf("Unable to load \"%s\" from %s\n", key.c_str(), filename.c_str());
    }
    return false;
  }
}

// Loads array from a .yaml file into an Eigen-type matrix or vector.
// Author: James Jackson
template <typename T>
bool get_yaml_eigen(const std::string key, const std::string filename, Eigen::MatrixBase<T>& val) 
{
  YAML::Node node = YAML::LoadFile(filename);
  std::vector<double> vec;
  if (node[key])
  {
    vec = node[key].as<std::vector<double>>();
    if (vec.size() == (val.rows() * val.cols()))
    {
      int k = 0;
      for (int i = 0; i < val.rows(); i++)
      {
        for (int j = 0; j < val.cols(); j++)
        {
          val(i,j) = vec[k++];
        }
      }
      return true;
    }
    else
    {
      printf("Eigen Matrix Size does not match parameter size for \"%s\" in %s", key.c_str(), filename.c_str());
      return false;
    }
  }
  else
  {
    printf("Unable to load \"%s\" from %s\n", key.c_str(), filename.c_str());
    return false;
  }
}

// Saturates a scalar value.
template <typename T>
T saturate(const T val, const T max, const T min)
{
  if (val > max)
    return max;
  if (val < min)
    return min;
  return val;
}

} // namespace common

#endif // COMMON_H
