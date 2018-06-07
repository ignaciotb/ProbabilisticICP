cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(icp_pkg)

find_package(PCL 1.5 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (icp_simple 
	src/icp_simple.cpp)

target_link_libraries (icp_simple ${PCL_LIBRARIES})