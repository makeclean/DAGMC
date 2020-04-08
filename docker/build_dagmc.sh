#!/bin/bash

# DAGMC source code already exists because it was copied into
# the docker image by travis.yml


set -ex

source ${docker_env}

if [ "${TRAVIS_PULL_REQUEST}" == "false" ]; then
  build_mw_reg_tests=ON
else
  build_mw_reg_tests=OFF
fi

rm -rf ${dagmc_build_dir}/bld ${dagmc_install_dir}
mkdir -p ${dagmc_build_dir}/bld
cd ${dagmc_build_dir}
cd bld
cmake ../DAGMC -DMOAB_DIR=${moab_install_dir} \
             -DBUILD_GEANT4=ON \
             -DGEANT4_DIR=${geant4_install_dir} \
             -DBUILD_CI_TESTS=ON \
             -DBUILD_MW_REG_TESTS=${build_mw_reg_tests} \
             -DBUILD_STATIC_EXE=${build_static_exe} \
             -DCMAKE_C_COMPILER=${CC} \
             -DCMAKE_CXX_COMPILER=${CXX} \
             -DCMAKE_Fortran_COMPILER=${FC} \
             -DCMAKE_INSTALL_PREFIX=${dagmc_install_dir}
make -j${jobs}
make install
cd
