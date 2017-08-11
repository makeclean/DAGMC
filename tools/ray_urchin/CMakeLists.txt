message(STATUS "Building ray_urchin")

# Source files
set(SRC_FILES)
list(APPEND SRC_FILES ray_urchin.cpp)

# Include directories
include_directories(.)
include_directories(${MOAB_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/dagmc)

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})
if (BUILD_STATIC)
   list(APPEND LINK_LIBS dagmc-static)
else()
   list(APPEND LINK_LIBS dagmc)
endif( BUILD_STATIC)

# Define executables
add_executable(ray_urchin  ${SRC_FILES})

# Link executables to moab/dagmc libraries
target_link_libraries(ray_urchin ${LINK_LIBS})
target_compile_features(ray_urchin PRIVATE cxx_range_for)
# Install
install(TARGETS ray_urchin DESTINATION ${INSTALL_BIN_DIR} PERMISSIONS ${PERMS})

