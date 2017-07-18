cmake_minimum_required(VERSION 3.0)
project(himalaya CXX Fortran)

#SET (CMAKE_Fortran_COMPILER gfortran)

#set the build type to release. If one needs debugging symbols, compile
#with "CMAKE_BUILD_TYPE Debug"
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set the CXX flags
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

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
add_library(himalaya STATIC ${LIBSOURCES})

#set the executable and specify the c++ version
add_executable(example ${CMAKE_CURRENT_SOURCE_DIR}/source/example.cpp)
target_compile_features(example PRIVATE cxx_range_for)

#link mh2l and himalaya to the executable
target_link_libraries(example himalaya)

#install
install(TARGETS example RUNTIME DESTINATION bin)
