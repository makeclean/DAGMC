#!/bin/bash
set -e 

# build hdf5
if [ ! -d $HOME/.ccache/hdf5 ] ; then
    wget https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.13/src/hdf5-1.8.13.tar.gz
    tar -zxf hdf5-1.8.13.tar.gz
    mv hdf5-1.8.13 hdf5
    cd hdf5
    mkdir bld
    cd bld
    ../configure --enable-shared --disable-debug --enable-optimize --prefix=$HOME/.ccache/hdf5
    make -j2 > log-file 2>&1
    make install
    cd ..
    cd ..
fi
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/.ccache/hdf5/lib/"
export PATH="$HOME/.ccache/hdf5/bin:$PATH"
