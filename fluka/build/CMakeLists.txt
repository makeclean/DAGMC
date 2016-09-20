message(STATUS "Building FluDAG")

enable_language(Fortran)
if(NASA)
  add_definitions(-DNASA)
  # if nasa build pull the GCR source information
  include(ExternalProject)
  ExternalProject_Add(sragcodes
  GIT_REPOSITORY https://github.com/kerrylee01/SRAGCodes
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  )
endif(NASA)

# Source files
set(SRC_FILES mainFluDAG.cpp)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS fludag)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(.)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# custom source routine
if(NASA)
  set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/bld/fluka/build/sragcodes-prefix/src/sragcodes/RadSource)
  set(FLUKA_FILE ${SOURCE_DIR}/source.f)
  set(FLUKA_INCLUDES ${FLUKA_DIR}/flukapro)

  # these are the hoops we have to jump through to link against the source routine
  add_custom_command(OUTPUT source.o
       COMMAND gfortran -c ${FLUKA_FILE} -I${FLUKA_INCLUDES} 
       DEPENDS ${FLUKA_FILE} sragcodes
       COMMENT "Custom command to compile source object with parenthetically named includes"
       )
#  add_custom_target(source DEPENDS source.o sraglib )
endif(NASA)

# Build FluDAG
if(NASA)
 add_executable (mainfludag ${SRC_FILES} source.o)
 target_link_libraries(mainfludag gcrsource ${LINK_LIBS})
else(NASA)
 add_executable (mainfludag ${SRC_FILES})
 target_link_libraries(mainfludag ${LINK_LIBS})
endif(NASA)


install(TARGETS mainfludag DESTINATION ${INSTALL_BIN_DIR}
                           PERMISSIONS ${PERMS})

