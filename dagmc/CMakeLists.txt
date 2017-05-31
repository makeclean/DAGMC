message(STATUS "Building DAGMC library")

# Source files
list(APPEND DAGMC_SRC_FILES dagmcmetadata.cpp DagMC.cpp)
                                
# Public headers
list(APPEND DAGMC_PUB_HEADERS dagmcmetadata.hpp DagMC.hpp)

# Include directories
include_directories(.)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Build dagmciface library
add_library(dagmc ${LIB_TYPE} ${DAGMC_SRC_FILES})
set_target_properties(dagmc PROPERTIES PUBLIC_HEADER "${DAGMC_PUB_HEADERS}")
install(TARGETS dagmc      LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                           ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                           PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
                           LIBRARY       PERMISSIONS ${PERMS}
                           ARCHIVE       PERMISSIONS ${PERMS})
list(APPEND LINK_LIBS pyne_dagmc)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS ${HDF5_LINK_LIBS})
target_link_libraries(dagmc ${LINK_LIBS})
add_subdirectory(test)
