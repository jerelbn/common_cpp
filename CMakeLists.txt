cmake_minimum_required(VERSION 3.5.1)
project(common_cpp)

# default to release build
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

find_package(Eigen3 REQUIRED)
find_package(GTest REQUIRED)

include_directories(include
  ${EIGEN3_INCLUDE_DIRS}
)

add_library(common_cpp INTERFACE)

add_executable(quaternion_test
  test/test_quaternion.cpp
)
target_link_libraries(quaternion_test
  common_cpp
  gtest
  gtest_main
  pthread
  ${GTEST_LIBRARIES}
)

add_executable(transform_test
  test/test_transform.cpp
)
target_link_libraries(transform_test
  common_cpp
  gtest
  gtest_main
  pthread
  ${GTEST_LIBRARIES}
)

add_executable(rotation_test
  test/test_rotation.cpp
)
target_link_libraries(rotation_test
  common_cpp
  gtest
  gtest_main
  pthread
  ${GTEST_LIBRARIES}
)

add_executable(matrix_test
  test/test_matrix.cpp
)
target_link_libraries(matrix_test
  common_cpp
  gtest
  gtest_main
  pthread
  ${GTEST_LIBRARIES}
)
