#!/bin/bash

# Make sure all the paths are correct

if [ -e build ];then
  rm -r build
fi

if [ -e install ];then
  rm -r install
fi

MY_PATH=$(pwd)

export CC=icx
export FC=ifx
export CXX=icpx

export BLAS_VENDOR=${BLAS_VENDOR:=MKL}
export BML_OPENMP=${BML_OPENMP:=yes}
export BML_GPU=${BML_GPU:=no}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export BML_TESTING=${BML_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export EXTRA_CFLAGS=${EXTRA_CFLAGS:=""}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:=""}
export CMAKE_Fortran_FLAGS="-fopenmp-targets=spir64 -DINTEL_SDK -fiopenmp"
export CMAKE_C_FLAGS="-fopenmp-targets=spir64 -DINTEL_SDK -D__STRICT_ANSI__ -DUSE_OMP_OFFLOAD -fiopenmp -Wno-missing-prototype-for-cc"

./build.sh configure

pushd build
make -j8
make install
popd
