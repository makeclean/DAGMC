message(STATUS "Building dagmc_povray")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES dagmc_povray.cpp)
list(APPEND SRC_FILES color_space.cpp)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Define executables
add_executable(dagmc_povray ${SRC_FILES})

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})

# Link executables to moab/dagmc libraries
target_link_libraries(dagmc_povray ${LINK_LIBS})

# Install
install(TARGETS dagmc_provray DESTINATION ${INSTALL_BIN_DIR} PERMISSIONS ${PERMS})

# build tests
#add_subdirectory(tests)