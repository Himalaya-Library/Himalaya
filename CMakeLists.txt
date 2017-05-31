cmake_minimum_required(VERSION 3.0)
project(h3cc CXX Fortran)

#SET (CMAKE_Fortran_COMPILER gfortran)

#Mac users should use g++ and modify the line below
#SET (CMAKE_CXX_COMPILER g++)

#include c++ header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source/include/hierarchies)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source/Eigen/)

#h3cc files
file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)

# mh2l files
file(GLOB_RECURSE FSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/include/mh2l/*f)

#build mh2l library and set linking language
ADD_LIBRARY(mh2l STATIC ${FSOURCES})
SET_TARGET_PROPERTIES(mh2l PROPERTIES LINKER_LANGUAGE Fortran)

#set the executable and specify the c++ version
add_executable(h3cc ${SOURCES})
target_compile_features(h3cc PRIVATE cxx_range_for)

#link mh2l to h3cc
TARGET_LINK_LIBRARIES(h3cc mh2l)

#install
install(TARGETS h3cc RUNTIME DESTINATION bin)

