cmake_minimum_required(VERSION 2.8.3)
project( common_cpp )

# choose debug or release build
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)

# enable multi-threading for Eigen
set(CMAKE_CXX_FLAGS "-std=c++0x -fopenmp -O3")

find_package( OpenCV REQUIRED )
find_package( Eigen3 REQUIRED )

include_directories(include)
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})
