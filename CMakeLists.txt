cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(dr_sasa)

FIND_PACKAGE( OpenMP )
  if(OPENMP_FOUND)
    #message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
#add_definitions(--std=c++11 -g)

file(GLOB_RECURSE files "src/*.cpp")

add_executable (dr_sasa ${files})
set_property(TARGET dr_sasa PROPERTY CXX_STANDARD 11)
set_property(TARGET dr_sasa PROPERTY CXX_STANDARD_REQUIRED ON)
add_definitions(-static -Wall -O3 -mcpu=native -mtune=native)

#target_link_libraries (dr_sasa ${OpenCL_LIBRARIES} -lgomp -lpthread)

#target_link_libraries (dr_sasa ${PCL_LIBRARIES} ${OpenCL_LIBRARIES})
