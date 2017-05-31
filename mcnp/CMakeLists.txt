message(STATUS "Including MCNP in dagmc library")

include_directories(${PROJECT_SOURCE_DIR}/dagmc)
include_directories(${PROJECT_SOURCE_DIR}/tally)
include_directories(${MOAB_INCLUDE_DIRS})

# dagmciface source files
set(MCNPFUNCS_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/meshtal_funcs.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/mcnp_funcs.cpp)

# dagmciface public headers
set(MCNPFUNCS_PUB_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/meshtal_funcs.h
			  ${CMAKE_CURRENT_SOURCE_DIR}/mcnp_funcs.h)

add_library(mcnp ${MCNPFUNCS_SRC_FILES})
target_link_libraries(mcnp ${DAG_LINK_LIBS} c++)
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "stdc++")

# Build MCNP5
if (BUILD_MCNP5)
  add_subdirectory(mcnp5)
endif (BUILD_MCNP5)

# Build MCNP6
if (BUILD_MCNP6)
  add_subdirectory(mcnp6)
endif (BUILD_MCNP6)

# Build tests
if (BUILD_MCNP5 OR BUILD_MCNP6)
  add_subdirectory(test)
endif (BUILD_MCNP5 OR BUILD_MCNP6)
