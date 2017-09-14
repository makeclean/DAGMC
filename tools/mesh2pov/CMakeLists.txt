message(STATUS "Building mesh2pov")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES mesh2pov.cpp)

# Libraries
set(LINK_LIBS ${MOAB_LIBRARIES})

# Define executables
add_executable(mesh2pov ${SRC_FILES})

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})

target_link_libraries(mesh2pov        ${LINK_LIBS})
install(TARGETS mesh2pov  DESTINATION ${INSTALL_BIN_DIR})

