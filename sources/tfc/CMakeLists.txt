# Cmake file for otis
cmake_minimum_required(VERSION 2.8.4)

# need c++11
add_definitions(-std=c++11)

IF(STATIC_LIB)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
ENDIF(STATIC_LIB)

# set the library
add_library(source_tf source_tf.cpp)

# and the includes
target_include_directories (source_tf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# global variable for cludgy Fluka building
SET_PROPERTY(GLOBAL PROPERTY TFC_FLUKA_SRC "${CMAKE_CURRENT_SOURCE_DIR}/source.f")
