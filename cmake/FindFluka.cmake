#find_path(FLUKA_LIBRARIES
#  NAMES libfluka.a libflukahp.a
#  HINTS ${FLUKA_DIR} ${FLUKA_DIR}/lib
#  PATHS ENV FLUKA_DIR
#  NO_DEFAULT_PATH
#  )

set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

find_library(FLUKA_LIBRARIES
  NAMES fluka flukahp
  PATHS ${FLUKA_DIR}
        ${FLUKA_DIR}/lib
  ENV ${FLUPRO}
)

if (FLUKA_LIBRARIES)
  set(FLUKA_LIBRARIES ${FLUKA_LIBRARIES} gfortran)
endif()

message(STATUS "FLUKA_LIBRARIES: ${FLUKA_LIBRARIES}")

if (FLUKA_LIBRARIES)
  message(STATUS "Found Fluka")
else ()
  message(FATAL_ERROR "Could not find Fluka")
endif ()
