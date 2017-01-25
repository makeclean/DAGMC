# Try to find Geant4
#
# Once done this will define
#
#  MOAB_FOUND - system has MOAB
#  MOAB_INCLUDE_DIRS - the MOAB include directory
#  MOAB_LIBRARIES - Link these to use MOAB
#  MOAB_DEFINITIONS - Compiler switches required for using MOAB

# GEANT4_FOUND - System has Geant4
# GEANT4_INCLUDE_DIR - The Geant4 include dir
# GEANT4_LIBRARIES - In case they are ever needed

find_path( GEANT4_CMAKE_CONFIG NAMES Geant4Config.cmake
  HINTS ${GEANT4_DIR}/lib/Geant4-10.0.2 ${GEANT4_DIR}/lib/Geant4-10.1.2 ${GEANT4_DIR}/lib/Geant4-10.2.1 ${GEANT4_DIR}
  PATHS ENV LD_LIBRARY_PATH DYLD_LIBRARY_PATH
  PATH_SUFFIXES lib Lib
  NO_DEFAULT_PATH
  )

MESSAGE ( STATUS "Found Geant4 in ${GEANT4_CMAKE_CONFIG}" )

INCLUDE ( ${GEANT4_CMAKE_CONFIG}/Geant4Config.cmake )

# Try to find Geant4Config.cmake
file(GLOB SEARCH_DIRS "${GEANT4_DIR}/lib*/Geant4-*")
string(REPLACE "\n" ";" SEARCH_DIRS ${SEARCH_DIRS})
find_path(GEANT4_CMAKE_CONFIG
          NAMES Geant4Config.cmake
          PATHS ${SEARCH_DIRS}
          NO_DEFAULT_PATH)
if (${GEANT4_CMAKE_CONFIG} MATCHES GEANT4_CMAKE_CONFIG-NOTFOUND)
  message(FATAL_ERROR "Geant4Config.cmake not found. Make sure to set -DGEANT4_DIR when running CMake.")
endif (${GEANT4_CMAKE_CONFIG} MATCHES GEANT4_CMAKE_CONFIG-NOTFOUND)
message(STATUS "Found Geant4Config.cmake at ${GEANT4_CMAKE_CONFIG}")

# Geant4 CMake config
include(${GEANT4_CMAKE_CONFIG}/Geant4Config.cmake)
