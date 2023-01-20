cmake_minimum_required(VERSION 3.0.2)
project(sfuise)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/" ${CMAKE_MODULE_PATH})

find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  nav_msgs
  sensor_msgs
  geometry_msgs
  message_generation
  cmake_modules
  cf_msgs
  isas_msgs
  sfuise_msgs
)

find_package(Eigen3 REQUIRED)
find_package(SuiteSparse REQUIRED MODULE)

catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS roscpp std_msgs message_runtime cf_msgs isas_msgs sfuise_msgs
  DEPENDS EIGEN3
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${BOOST_INCLUDE_DIR}
  ${SUITESPARSE_INCLUDE_DIRS}
)

add_executable(EstimationInterface src/EstimationInterface.cpp)
target_link_libraries(EstimationInterface ${catkin_LIBRARIES})

add_executable(SplineFusion src/SplineFusion.cpp)
target_link_libraries(SplineFusion ${catkin_LIBRARIES} cholmod)

set(CMAKE_BUILD_TYPE "Debug")
