cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(dr_sasa)
# find_package( PCL 1.2 REQUIRED )
# include_directories(${PCL_INCLUDE_DIRS})
# link_directories(${PCL_LIBRARY_DIRS})
# add_definitions(${PCL_DEFINITIONS})

#find_package( OpenCL REQUIRED )
#include_directories(${OPENCL_INCLUDE_DIRS})

find_package(OpenMP REQUIRED)

add_definitions(--std=c++11 -static ${OpenMP_CXX_FLAGS} -Wall -Ofast -mcpu=native)
#add_definitions(--std=c++11 -g)

file(GLOB_RECURSE files "src/*.cpp")

add_executable (dr_sasa ${files})
target_link_libraries (dr_sasa -lgomp -lpthread)
#target_link_libraries (dr_sasa ${OpenCL_LIBRARIES} -lgomp -lpthread)

#target_link_libraries (dr_sasa ${PCL_LIBRARIES} ${OpenCL_LIBRARIES})
