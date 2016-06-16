#!/bin/bash
set -e 

# build moab
if [ ! -d $HOME/.ccache/moab ] ; then
    git clone https://bitbucket.org/fathomteam/moab
    cd moab
    git checkout $MOAB_VERSION
    autoreconf -fi
    mkdir bld
    cd bld
    ../configure --enable-dagmc --enable-shared --disable-debug --enable-optimize --with-hdf5=$HOME/.ccache/hdf5 --prefix=$HOME/.ccache/moab
    make -j2 > log-file 2>&1
    make install
    cd ..
    cd ..
fi
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/.ccache/moab/lib/"
export PATH="$HOME/.ccache/moab/bin:$PATH"
