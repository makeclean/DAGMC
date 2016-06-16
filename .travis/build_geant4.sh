#!/bin/bash
set -e 

# build geant4
if [ ! -d $HOME/.ccache/geant4.10.00.p02 ] ; then
    wget http://geant4.cern.ch/support/source/geant4.10.00.p02.tar.gz
    tar -zxf geant4.10.00.p02.tar.gz
    cd geant4.10.00.p02
    mkdir bld
    cd bld
    cmake ../. -DCMAKE_INSTALL_PREFIX=$HOME/.ccache/geant4.10.00.p02
    # geant4 compile takes 30 mins so need to use travis_wait
    travis_wait make -j2
    make install
    cd ..
    cd ..
fi
export PATH="$HOME/.ccache/geant4.10.00.p02/bin:$PATH"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/.ccache/geant4.10.00.p02/lib/"
