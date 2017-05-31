message(STATUS "Including UWUW in dagmc library")

# dagmc source files
set(DAGMC_SRC_FILES ${DAGMC_SRC_FILES} ${CMAKE_CURRENT_LIST_DIR}/uwuw.cpp PARENT_SCOPE)

# dagmc public headers
set(DAGMC_PUB_HEADERS ${DAGMC_PUB_HEADERS} ${CMAKE_CURRENT_LIST_DIR}/uwuw.hpp PARENT_SCOPE)

# Include directories
include_directories(..)
include_directories(${CMAKE_SOURCE_DIR}/dagmc)
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MOAB_INCLUDE_DIRS})

set (UWUW_PREPROC_SRC "uwuw_preproc.cpp")
add_executable(uwuw_preproc ${UWUW_PREPROC_SRC})
target_link_libraries(uwuw_preproc dagmc  
  ${DAG_LINK_LIBS} ${HDF5_LINK_LIBS}  
  uwuw_preprocessor pyne_dagmc)

set(UWUW_PREPROCESSOR_SRC "uwuw_preprocessor.cpp")
add_library(uwuw_preprocessor ${UWUW_PREPROCESSOR_SRC})
target_link_libraries(uwuw_preprocessor dagmc
  ${DAG_LINK_LIBS} ${HDF5_LINK_LIBS}
  pyne_dagmc)

# Build UWUW tests
add_subdirectory(tests)

install(TARGETS uwuw_preproc DESTINATION ${INSTALL_BIN_DIR} PERMISSIONS ${PERMS})
