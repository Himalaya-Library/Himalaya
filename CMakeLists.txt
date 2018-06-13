cmake_minimum_required(VERSION 3.1)
project(Himalaya CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# require an Eigen installation
find_package(Eigen3 3.0 REQUIRED)
find_package(Mathematica 8.0)

# if needed, set your Fortran and C++ compiler
#set (CMAKE_Fortran_COMPILER gfortran)
#set (CMAKE_CXX_COMPILER clang++-3.8)
#set (CMAKE_CXX_COMPILER icpc)

# set the build type to release. If one needs debugging symbols, compile
# with "CMAKE_BUILD_TYPE Debug"
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set version number
set(Himalaya_VERSION_MAJOR 2)
set(Himalaya_VERSION_MINOR 0)
set(Himalaya_VERSION_RELEASE 0)

# set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set the CXX flags
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

# set paths
set(INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/source/include)

#include c++ header files
include_directories(${INCLUDE_PATH})
include_directories(${INCLUDE_PATH}/hierarchies)
include_directories(${INCLUDE_PATH}/mh2_eft)
include_directories(${INCLUDE_PATH}/misc)
include_directories(${EIGEN3_INCLUDE_DIR})

# setting the correct version to all headers
configure_file (
  "${INCLUDE_PATH}/version.hpp.in"
  "${INCLUDE_PATH}/version.hpp"
  )

# Himalaya library files
file(GLOB SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/source/hierarchies/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/source/mh2_eft/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/source/misc/*.cpp)
set(LIBSOURCES ${SOURCES})

# DSZ library files
file(GLOB FSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/include/mh2l/*.f)

# build the himalaya and DSZ library
file(GLOB_RECURSE HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.hpp)
file(GLOB LIBSOURCES ${SOURCES} ${HEADERS})
list(REMOVE_ITEM LIBSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/example.cpp)
add_library(DSZ_static STATIC ${FSOURCES})
set_target_properties(DSZ_static PROPERTIES OUTPUT_NAME DSZ)
add_library(DSZ_shared SHARED ${FSOURCES})
set_target_properties(DSZ_shared PROPERTIES OUTPUT_NAME DSZ)
add_library(Himalaya_static STATIC ${LIBSOURCES})
set_target_properties(Himalaya_static PROPERTIES OUTPUT_NAME Himalaya)
add_library(Himalaya_shared SHARED ${LIBSOURCES})
set_target_properties(Himalaya_shared PROPERTIES OUTPUT_NAME Himalaya)

# set the executable
add_executable(example ${CMAKE_CURRENT_SOURCE_DIR}/source/example.cpp)

# link mh2l and himalaya to the executable
target_link_libraries(example Himalaya_static DSZ_static)

# create LibrayLink executable
if(Mathematica_FOUND)
  set(LL_SRC ${CMAKE_CURRENT_SOURCE_DIR}/source/LibraryLink/Himalaya_LibraryLink.cpp)
  set(LL_LIB Himalaya_LibraryLink)
  include_directories(${Mathematica_INCLUDE_DIR})
  include_directories(${Mathematica_MathLink_INCLUDE_DIR})
  Mathematica_ADD_LIBRARY(${LL_LIB} ${LL_SRC})
  target_link_libraries(${LL_LIB} PRIVATE Himalaya_shared DSZ_shared ${Mathematica_MathLink_LIBRARIES})
  set_target_properties(${LL_LIB} PROPERTIES LINK_FLAGS "${Mathematica_MathLink_LINKER_FLAGS}")
  Mathematica_ABSOLUTIZE_LIBRARY_DEPENDENCIES(${LL_LIB})
endif()

# install
install(
  TARGETS Himalaya_static Himalaya_shared
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)

# adding documentation
add_subdirectory(doc)

# test
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test)
  enable_testing()
  add_subdirectory(test)
endif()
