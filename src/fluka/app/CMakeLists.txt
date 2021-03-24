enable_language(Fortran)
# add any files that you would like to use to build a custom fluka
# e.g. add source.f and so on
set(FLUKA_FILES bfield.cpp magfld.f plasma_source.cpp source.f)
set(SRC_FILES mainFluDAG.cpp ${FLUKA_FILES})

set(LINK_LIBS fludag)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/src/fluka)
include_directories(${FLUKA_DIR}/include)

dagmc_install_exe(mainfludag)
