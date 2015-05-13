# Cmake file for otis
cmake_minimum_required(VERSION 2.8.4)

# need c++11
add_definitions(-std=c++11)

# set the library
add_library(source_tf source_tf.cpp)

# and the includes
target_include_directories (source_tf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# global variable for cludgy Fluka building
SET_PROPERTY(GLOBAL PROPERTY TFC_FLUKA_SRC "${CMAKE_CURRENT_SOURCE_DIR}/source.f")
