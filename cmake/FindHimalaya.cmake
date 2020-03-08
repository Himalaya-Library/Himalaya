# Distributed under the GPLv3 License.
# Author: Alexander Voigt
#
# FindHimalaya
# -------------
#
# Finds the Himalaya library.
# [https://github.com/Himalaya-Library/Himalaya]
#
# This module reads the following variables:
#
# Himalaya_LIBRARY       - Himalaya library directory
# Himalaya_INCLUDE_DIR   - Himalaya include directory
#
# This module defines the following variables:
#
# Himalaya_FOUND         - set if Himalaya has been found
# Himalaya_VERSION       - Himalaya version
# Himalaya_INCLUDE_DIRS  - Himalaya include directory
# Himalaya_LIBRARIES     - Himalaya library
# DSZ_LIBRARIES          - Library with MSSM 2-loop corrections O(at*as)
#                          G. Degrassi, P. Slavich and F. Zwirner, 
#                          Nucl. Phys. B611 (2001) 403 [hep-ph/0105096]
#
# and defines the following imported targets:
#
# Himalaya::Himalaya
# Himalaya::DSZ

# search himalaya/Himalaya_interface.hpp first in ${Himalaya_INCLUDE_DIR}
find_path(Himalaya_INCLUDE_DIRS
  NAMES himalaya/Himalaya_interface.hpp
  PATHS
    ${Himalaya_INCLUDE_DIR}
  PATH_SUFFIXES
    include
  NO_DEFAULT_PATH
  NO_CMAKE_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  NO_CMAKE_FIND_ROOT_PATH
)

find_path(Himalaya_INCLUDE_DIRS
  NAMES himalaya/Himalaya_interface.hpp
)

# search Himalaya library first in ${Himalaya_LIBRARY}
find_library(Himalaya_LIBRARIES
  NAMES Himalaya
  PATHS
    ${Himalaya_LIBRARY}
  PATH_SUFFIXES
    build
  NO_DEFAULT_PATH
  NO_CMAKE_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  NO_CMAKE_FIND_ROOT_PATH
)

find_library(Himalaya_LIBRARIES
  NAMES Himalaya
)

# search DSZ library first in ${Himalaya_LIBRARY}
find_library(DSZ_LIBRARIES
  NAMES DSZ
  PATHS
    ${Himalaya_LIBRARY}
  PATH_SUFFIXES
    build
  NO_DEFAULT_PATH
  NO_CMAKE_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  NO_CMAKE_FIND_ROOT_PATH
)

find_library(DSZ_LIBRARIES
  NAMES DSZ
)

# find version
if(Himalaya_INCLUDE_DIRS)
  file(READ "${Himalaya_INCLUDE_DIRS}/himalaya/version.hpp" _himalaya_version_header)

  string(REGEX MATCH "define[ \t]+Himalaya_VERSION_MAJOR[ \t]+([0-9]+)" _himalaya_version_major_match "${_himalaya_version_header}")
  set(Himalaya_VERSION_MAJOR "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+Himalaya_VERSION_MINOR[ \t]+([0-9]+)" _himalaya_version_minor_match "${_himalaya_version_header}")
  set(Himalaya_VERSION_MINOR "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+Himalaya_VERSION_RELEASE[ \t]+([0-9]+)" _himalaya_version_release_match "${_himalaya_version_header}")
  set(Himalaya_VERSION_RELEASE "${CMAKE_MATCH_1}")

  set(Himalaya_VERSION ${Himalaya_VERSION_MAJOR}.${Himalaya_VERSION_MINOR}.${Himalaya_VERSION_RELEASE})

  if(Himalaya_FIND_VERSION)
    if(Himalaya_FIND_VERSION_EXACT AND NOT ${Himalaya_VERSION} VERSION_EQUAL ${Himalaya_FIND_VERSION})
      message(FATAL_ERROR "Himalaya version ${Himalaya_VERSION} found in ${Himalaya_INCLUDE_DIRS}, "
        "but exact version ${Himalaya_FIND_VERSION} is required.")
    elseif(${Himalaya_VERSION} VERSION_LESS ${Himalaya_FIND_VERSION})
      message(FATAL_ERROR "Himalaya version ${Himalaya_VERSION} found in ${Himalaya_INCLUDE_DIRS}, "
        "but at least version ${Himalaya_FIND_VERSION} is required.")
    endif()
  endif()
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Himalaya
  FOUND_VAR Himalaya_FOUND
  REQUIRED_VARS
    Himalaya_LIBRARIES
    Himalaya_INCLUDE_DIRS
    DSZ_LIBRARIES
)

if(Himalaya_FOUND AND NOT TARGET Himalaya::Himalaya)
  add_library(Himalaya::Himalaya UNKNOWN IMPORTED)
  set_target_properties(Himalaya::Himalaya PROPERTIES
    IMPORTED_LOCATION "${Himalaya_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${Himalaya_INCLUDE_DIRS}"
  )
endif()

if(Himalaya_FOUND AND NOT TARGET Himalaya::DSZ)
  add_library(Himalaya::DSZ UNKNOWN IMPORTED)
  set_target_properties(Himalaya::DSZ PROPERTIES
    IMPORTED_LOCATION "${DSZ_LIBRARIES}"
  )
endif()

mark_as_advanced(
  Himalaya_INCLUDE_DIRS
  Himalaya_LIBRARIES
  DSZ_LIBRARIES
)
