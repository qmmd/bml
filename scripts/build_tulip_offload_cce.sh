#!/bin/bash

# Make sure all the paths are correct

rm -r build
rm -r install

MY_PATH=$(pwd)

export CC=${CC:=cc}
export FC=${FC:=ftn}
export CXX=${CXX:=CC}
export BLAS_VENDOR=${BLAS_VENDOR:=Auto}
export BML_OPENMP=${BML_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export BML_MAGMA=${BML_MAGMA:=yes}
export MAGMA_ROOT=${MAGMA_HOME}
export BML_TESTING=${BML_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS:="-ef -DCRAY_SDK"}
#export EXTRA_CFLAGS=${EXTRA_CFLAGS:="-Ofast -DUSE_OMP_OFFLOAD -DCRAY_SDK"}
export CMAKE_C_FLAGS=${CMAKE_C_FLAGS:="-Ofast -DUSE_OMP_OFFLOAD -DCRAY_SDK"}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-L${BLASDIR} -l${BLASLIB}"}

./build.sh configure

pushd build
make -j8
make install
popd
