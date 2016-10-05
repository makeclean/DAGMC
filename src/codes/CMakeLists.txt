# This CMakeFile determines which of the physics code
# interfaces gets built

if (BUILD_GEANT4)
   add_subdirectory(geant4)
endif (BUILD_GEANT4)

if (BUILD_FLUKA)
   add_subdirectory(fluka)
endif (BUILD_FLUKA)

if (BUILD_MCNP5) 
   add_subdirectory(mcnp)
endif(BUILD_MCNP5)

