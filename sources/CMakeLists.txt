message(STATUS "Building Source Sampling library")

# Source files
file(GLOB SRC_FILES ${SRC_FILES} "*.cpp")

# Public headers
file(GLOB PUB_HEADERS ${PUB_HEADERS} "*.h")

# Linked libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
#list(APPEND LINK_LIBS ${DF5_hdf5_hl_LIBRARY})
#list(APPEND LINK_LIBS ${HDF5_hdf5_LIBRARY})

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})
#include_directories(${HDF5_INCLUDE_DIRS})

# Build pyne_dagmc library
add_library(source_sampling ${LIB_TYPE} ${SRC_FILES})
#set_target_properties(source_sampling PROPERTIES PUBLIC_HEADER PUB_HEADERS)
target_link_libraries(source_sampling ${LINK_LIBS})
#install(TARGETS source_sampling LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
#                           ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
#                           PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
#                           LIBRARY       PERMISSIONS ${PERMS}
#                           ARCHIVE       PERMISSIONS ${PERMS})
