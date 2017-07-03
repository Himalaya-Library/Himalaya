cmake_minimum_required(VERSION 3.0)
project(himalaya CXX Fortran)

#SET (CMAKE_Fortran_COMPILER gfortran)

#set the CXX flags
SET(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -pedantic")

#include c++ header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source/include/hierarchies)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source/Eigen/)

#himalaya files
file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
set(LIBSOURCES ${SOURCES})

# mh2l files
file(GLOB_RECURSE FSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/include/mh2l/*f)

#build the himalaya library and combine it wtih mh2l
file(GLOB_RECURSE HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.hpp)
file(GLOB LIBSOURCES ${SOURCES} ${HEADERS} ${FSOURCES})
list(REMOVE_ITEM LIBSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/example.cpp)
ADD_LIBRARY(himalaya STATIC ${LIBSOURCES})

#set the executable and specify the c++ version
add_executable(example ${CMAKE_CURRENT_SOURCE_DIR}/source/example.cpp)
target_compile_features(example PRIVATE cxx_range_for)

#link mh2l and himalaya to the executable
TARGET_LINK_LIBRARIES(example himalaya)

#install
install(TARGETS example RUNTIME DESTINATION bin)
