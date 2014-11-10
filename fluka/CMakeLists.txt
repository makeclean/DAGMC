CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
MESSAGE ( STATUS "Including FLUKA in dagmciface library.")

SET ( CODE_SRCDIR      "fluka" )

SET ( CODE_SRC_FILES   "${CODE_SRCDIR}/fluka_funcs.cpp" )
SET ( CODE_PUB_HEADERS "${CODE_SRCDIR}/fluka_funcs.h" )

#SET ( CODE_SRC_FILES   /mnt/data/opt/mcnp5v16_current/Source/src/plasma_source.cpp )


SET ( CODE_CONFIG       )

SET ( DI_PUB_HEADERS ${DI_PUB_HEADERS}
                     ${CODE_PUB_HEADERS}
                     ${CODE_CONFIG}
      PARENT_SCOPE)

SET ( DI_SRC_FILES ${DI_SRC_FILES} 
                   ${CODE_SRC_FILES} 
      PARENT_SCOPE)

#INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/${CODE_SRCDIR})

INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${CODE_SRC_DIR} )

ENABLE_LANGUAGE( Fortran )

SET(FLUDAG_LIBS
     dagmciface
     ${HDF5_LIBRARIES}
     ${PYNE_LIBS_DIR}/libpyne.so
     ${DAGMC_LIBRARIES} 
     gfortran
  )

# If making mainfludag, define FLUPRO)
# add_custom_command (TARGET mainfludag
# PRE_BUILD COMMAND 

# custom sources
SET(FLUKA_FILE ${FLUKA_DIR}/usermvax/source.f)
SET (FLUKA_INCLUDES ${FLUKA_DIR}/flukapro)
ADD_CUSTOM_COMMAND(OUTPUT source.o
       COMMAND gfortran -c ${FLUKA_FILE} -I${FLUKA_INCLUDES} 
       DEPENDS ${FLUKA_DIR}/usermvax/source.f 
       COMMENT "Custom command to compile source object with parenthetically named includes"
       )

ADD_CUSTOM_TARGET (source DEPENDS source.o )
#/mnt/data/opt/mcnp5v16_current/Source/src/plasma_source.cpp
ADD_EXECUTABLE ( mainfludag mainFluDAG.cpp /mnt/data/opt/mcnp5v16_current/Source/src/plasma_source.cpp source.o)
INSTALL( TARGETS mainfludag 
  RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
# Order very important!
TARGET_LINK_LIBRARIES ( mainfludag ${FLUDAG_LIBS} ${FLUKA_DIR}/libflukahp.a gfortran)

ADD_SUBDIRECTORY(tests)

