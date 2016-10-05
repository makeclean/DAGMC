message(STATUS "Building dagmciface library")


# Source files
set(DI_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/pyne/pyne.cpp 
		 ${CMAKE_CURRENT_SOURCE_DIR}/uwuw/uwuw.cpp 
		 )

# Public headers
set(DI_PUB_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/pyne/pyne.h 
  		   ${CMAKE_CURRENT_SOURCE_DIR}/uwuw/uwuw.hpp 
		   )

# Include directories
include_directories(.)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

# Build dagmciface library
add_library(dagmciface ${LIB_TYPE} ${DI_SRC_FILES})

# Linked libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${HDF5_hdf5_hl_LIBRARY})
list(APPEND LINK_LIBS ${HDF5_hdf5_LIBRARY})

target_link_libraries(dagmciface ${LINK_LIBS})

# Include directories
include_directories(.)

set_target_properties(dagmciface PROPERTIES PUBLIC_HEADER "${DI_PUB_HEADERS}")
install(TARGETS dagmciface LIBRARY       DESTINATION ${INSTALL_LIB_DIR}
                           ARCHIVE       DESTINATION ${INSTALL_LIB_DIR}
                           PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
                           LIBRARY       PERMISSIONS ${PERMS}
                           ARCHIVE       PERMISSIONS ${PERMS})

add_subdirectory(tally)
