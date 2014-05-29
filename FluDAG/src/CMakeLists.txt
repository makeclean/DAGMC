PROJECT (mainFluDAG)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)
# CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

IF (NOT DEFINED HOME) 
	set (HOME "$ENV{HOME}")
ENDIF (NOT DEFINED HOME) 

IF (NOT DEFINED FLUPRO)
	set (FLUPRO "$ENV{FLUPRO}")
ENDIF (NOT DEFINED FLUPRO)

#############################################
# PROJECT (PyneDag)
#

SET(pynelibs
   ${HDF5_HOME}/lib/libhdf5.so
   ${PYNE_HOME}/lib/libpyne.so
   ${PYNE_HOME}/lib/libpyne_material.so
   ${PYNE_HOME}/lib/libpyne_nucname.so
   ${PYNE_HOME}/lib/libpyne_data.so
)

SET(moablibs 
	${MOAB_HOME}/lib/libdagmc.so 
	${MOAB_HOME}/lib/libMOAB.so
)

# Where to look for includes
INCLUDE_DIRECTORIES(${MOAB_HOME}/include ${HDF5_HOME}/include ${PYNE_HOME}/include)
ENABLE_LANGUAGE(Fortran)

#SET (fluka_file ${FLUPRO}/usermvax/mgdraw.f)
SET (fluka_includes ${FLUPRO}/flukapro)
SET (fluka_libdir ${FLUPRO})
#ADD_CUSTOM_COMMAND(OUTPUT mgdraw.o
#       COMMAND gfortran -c ${fluka_file} -I${fluka_includes} 
#       DEPENDS ${FLUPRO}/usermvax/mgdraw.f
#       COMMENT "Custom command to compile mgdraw object with parenthetically named includes"
#       )
#ADD_CUSTOM_TARGET (mgdraw
#	        DEPENDS mgdraw.o)

# Create the executable that can call flukam_
ADD_EXECUTABLE(mainfludag mainFluDAG.cpp fluka_funcs.cpp)
TARGET_LINK_LIBRARIES(mainfludag ${moablibs} ${pynelibs} ${FLUPRO}/libflukahp.a gfortran)
