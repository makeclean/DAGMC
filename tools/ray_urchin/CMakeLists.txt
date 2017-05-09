message(STATUS "Building ray_urchin")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES ray_urchin.cpp)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
list(APPEND LINK_LIBS dagmciface)

# Define executables
add_executable(ray_urchin  ${SRC_FILES})

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/common)

# Link executables to moab/dagmc libraries
target_link_libraries(ray_urchin ${LINK_LIBS})

# Install
install(TARGETS ray_urchin DESTINATION ${INSTALL_BIN_DIR} PERMISSIONS ${PERMS})

