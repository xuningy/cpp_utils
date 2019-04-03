cmake_minimum_required(VERSION 2.8.3)

project(cpp_utils)

find_package(catkin REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS
)

include_directories(include ${catkin_INCLUDE_DIRS})
link_directories(${catkin_LIBRARY_DIRS})
add_definitions(${catkin_DEFINITIONS})
