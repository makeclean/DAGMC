message(STATUS "Building DagSolid unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS gtest)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS ${MOAB_LIBRARIES})
list(APPEND LINK_LIBS dagmc)

# Include directories
include_directories(.)
include_directories(dagmc)
include_directories(${GTEST_INCLUDE_DIR})

# Setup unit tests
set(DRIVER dagmc_test_driver.cpp)
setup_test(dagmc_simple_test cpp ${DRIVER} "${LINK_LIBS}")

# Install h5m files
install(FILES test_geom.h5m DESTINATION tests)

# Enable unit tests
enable_testing()
