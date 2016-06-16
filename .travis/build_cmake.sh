#!/bin/bash
set -e 

# check to see if its a cached build
if [ ! -d $HOME/.ccache/cmake ] ; then
    git clone git://cmake.org/cmake.git custom_cmake
    cd custom_cmake
    git checkout master
    mkdir bld
    cd bld
    cmake ../. -DCMAKE_INSTALL_PREFIX=$HOME/.ccache/cmake
    make -j2 > log-file 2>&1
    make install
    cd ..
    cd .. 
fi
# always append to the path
export PATH="$HOME/.ccache/cmake/bin:$PATH"


