message(STATUS "Building DAGMC Library")

# Include directories
include_directories(.)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Build dagmc library
add_library(dagmc ${LIB_TYPE} DagMC.cpp)

# Linked libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${HDF5_hdf5_hl_LIBRARY})
list(APPEND LINK_LIBS ${HDF5_hdf5_LIBRARY})
list(APPEND LINK_LIBS ${MOAB_LIBRARIES})
# Include directories
include_directories(.)

target_link_libraries(dagmc ${LINK_LIBS})

set_target_properties(dagmc PROPERTIES PUBLIC_HEADER DagMC.hpp)
install(TARGETS dagmc      LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                           ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                           PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
                           LIBRARY       PERMISSIONS ${PERMS}
                           ARCHIVE       PERMISSIONS ${PERMS})
