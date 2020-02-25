cmake_minimum_required(VERSION 3.1)
project(Himalaya LANGUAGES CXX Fortran)

include(GNUInstallDirs)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# require an Eigen installation
find_package(Eigen3 3.0 REQUIRED)
find_package(Mathematica 8.0)

# if needed, set your Fortran and C++ compiler
#set (CMAKE_Fortran_COMPILER gfortran)
#set (CMAKE_CXX_COMPILER clang++)
#set (CMAKE_CXX_COMPILER icpc)

# set the build type to release. If one needs debugging symbols, compile
# with "CMAKE_BUILD_TYPE Debug"
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# enable all warnings in Debug mode
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  if(CMAKE_COMPILER_IS_GNUCXX OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
  endif()
endif()

# set paths
set(SOURCE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(INCLUDE_PATH ${SOURCE_PATH}/include)

#include c++ header files
include_directories(${INCLUDE_PATH}/hierarchies)
include_directories(${INCLUDE_PATH}/mh2_eft)
include_directories(${INCLUDE_PATH}/mh2_fo)
include_directories(${INCLUDE_PATH}/mh2l)
include_directories(${INCLUDE_PATH}/misc)

# generate himalaya.pc file for pkg-config
configure_file (
  "${PROJECT_SOURCE_DIR}/cmake/himalaya.pc.in"
  "${CMAKE_BINARY_DIR}/himalaya.pc"
  @ONLY)

# Himalaya library files
set(HIM_LIBSOURCES
  # interface
  ${SOURCE_PATH}/HierarchyCalculator.cpp
  ${SOURCE_PATH}/HierarchyObject.cpp
  ${SOURCE_PATH}/Himalaya_interface.cpp
  # hierarchies
  ${SOURCE_PATH}/hierarchies/H32q2g.cpp
  ${SOURCE_PATH}/hierarchies/H3.cpp
  ${SOURCE_PATH}/hierarchies/H3q22g.cpp
  ${SOURCE_PATH}/hierarchies/H4.cpp
  ${SOURCE_PATH}/hierarchies/H5.cpp
  ${SOURCE_PATH}/hierarchies/H5g1.cpp
  ${SOURCE_PATH}/hierarchies/H6b2qg2.cpp
  ${SOURCE_PATH}/hierarchies/H6b.cpp
  ${SOURCE_PATH}/hierarchies/H6bq22g.cpp
  ${SOURCE_PATH}/hierarchies/H6bq2g2.cpp
  ${SOURCE_PATH}/hierarchies/H6.cpp
  ${SOURCE_PATH}/hierarchies/H6g2.cpp
  ${SOURCE_PATH}/hierarchies/H9.cpp
  ${SOURCE_PATH}/hierarchies/H9q2.cpp
  # EFT calculation
  ${SOURCE_PATH}/mh2_eft/Mh2EFTCalculator.cpp
  ${SOURCE_PATH}/mh2_eft/ThresholdCalculator.cpp
  ${SOURCE_PATH}/mh2_eft/threshold_loop_functions.cpp
  # fixed-order calculation
  ${SOURCE_PATH}/mh2_fo/linalg2.cpp
  ${SOURCE_PATH}/mh2_fo/MSSM_mass_eigenstates.cpp
  ${SOURCE_PATH}/mh2_fo/pv.cpp
  # helpers
  ${SOURCE_PATH}/misc/dilog.cpp)

# DSZ library files
set(DSZ_LIBSOURCES
  ${SOURCE_PATH}/mh2l/DSZHiggs.cpp
  ${SOURCE_PATH}/mh2l/DSZHiggs.f
  ${SOURCE_PATH}/mh2l/functs.f
  ${SOURCE_PATH}/mh2l/li2.cpp)

# build the himalaya and DSZ library
add_library(DSZ_static STATIC ${DSZ_LIBSOURCES})
set_target_properties(DSZ_static PROPERTIES OUTPUT_NAME DSZ)
target_include_directories(DSZ_static
  PUBLIC
    ${EIGEN3_INCLUDE_DIR}
)

add_library(DSZ_shared SHARED ${DSZ_LIBSOURCES})
set_target_properties(DSZ_shared PROPERTIES OUTPUT_NAME DSZ)
target_include_directories(DSZ_shared
  PUBLIC
    ${EIGEN3_INCLUDE_DIR}
)

add_library(Himalaya_static STATIC ${HIM_LIBSOURCES})
set_target_properties(Himalaya_static PROPERTIES OUTPUT_NAME Himalaya)
target_include_directories(Himalaya_static
  PUBLIC
    ${INCLUDE_PATH}
    ${EIGEN3_INCLUDE_DIR}
)

add_library(Himalaya_shared SHARED ${HIM_LIBSOURCES})
set_target_properties(Himalaya_shared PROPERTIES OUTPUT_NAME Himalaya)
target_link_libraries(Himalaya_shared PUBLIC DSZ_shared)
target_include_directories(Himalaya_shared
  PUBLIC
    ${INCLUDE_PATH}
    ${EIGEN3_INCLUDE_DIR}
)

# set the executable
add_executable(example ${SOURCE_PATH}/example.cpp)
target_link_libraries(example Himalaya_static DSZ_static)

# create LibrayLink executable
if(Mathematica_FOUND)
  set(LL_SRC ${SOURCE_PATH}/LibraryLink/Himalaya_LibraryLink.cpp)
  set(LL_LIB Himalaya_LibraryLink)
  Mathematica_ADD_LIBRARY(${LL_LIB} ${LL_SRC})
  target_link_libraries(${LL_LIB} PRIVATE Himalaya_shared ${Mathematica_MathLink_LIBRARIES})
  target_include_directories(${LL_LIB}
    PRIVATE
      ${Mathematica_INCLUDE_DIR}
      ${Mathematica_MathLink_INCLUDE_DIR}
    )
  set_target_properties(${LL_LIB} PROPERTIES LINK_FLAGS "${Mathematica_MathLink_LINKER_FLAGS}")

  Mathematica_ABSOLUTIZE_LIBRARY_DEPENDENCIES(${LL_LIB})
endif()

# install libraries
install(
  TARGETS Himalaya_static Himalaya_shared DSZ_static DSZ_shared
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# install public headers
install(FILES
  "${INCLUDE_PATH}/HierarchyCalculator.hpp"
  "${INCLUDE_PATH}/HierarchyObject.hpp"
  "${INCLUDE_PATH}/Himalaya_interface.hpp"
  "${INCLUDE_PATH}/version.hpp"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/himalaya")

# install himalaya.pc file
install(FILES
  "${CMAKE_BINARY_DIR}/himalaya.pc"
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig")

# adding documentation
add_subdirectory(doc)

# test
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test)
  enable_testing()
  add_subdirectory(test)
endif()
