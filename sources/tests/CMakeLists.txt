# CMAKE build script for DAGMC Tally unit testing

# CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Where to look for includes
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/source )

MESSAGE(STATUS "${GTEST_HOME}")

# if we want a static executable
IF(STATIC_LIB)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
ENDIF(STATIC_LIB)

# set libraries
SET(LIBRARIES
    pthread
    plasma_sources
    gtest
    pthread
)


# set up DAGMC Tally test cases
ADD_EXECUTABLE(test_PlasmaSources test_PlasmaSources.cpp sources_unit_test_driver.cc)
TARGET_LINK_LIBRARIES(test_PlasmaSources ${LIBRARIES})

INSTALL(TARGETS test_PlasmaSources DESTINATION tests)

# files needed for tests
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/rzsource.txt" DESTINATION tests)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/r_prof.txt" DESTINATION tests)

# enable DAGMC Tally test cases
ENABLE_TESTING()

# list DAGMC Tally test cases
ADD_TEST(test_PlasmaSources test_PlasmaSources)
