cmake_minimum_required(VERSION 2.8.3)
project(common_cpp)

# default to release build
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
set(CMAKE_CXX_FLAGS "-std=c++11 -openmp")

## Find libraries
find_package(Eigen3 REQUIRED)
#find_package(OpenCV 3 REQUIRED)
#message(STATUS "OpenCV Version: ${OpenCV_VERSION} (${OpenCV_DIR})")

## Specify include directories
include_directories(include
  ${Eigen3_INCLUDE_DIRS}
  # ${OpenCV_INCLUDE_DIRS}
) 
