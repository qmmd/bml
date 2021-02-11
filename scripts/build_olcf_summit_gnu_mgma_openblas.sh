#!/bin/bash
module load cmake
module load cuda
module load gcc/8.1.1
module load netlib-lapack
module load openblas
module load magma

rm -r build
rm -r install

MY_PATH=$(pwd)

#get jsrun with full path
JSRUN=$(which jsrun)
echo ${JSRUN}

export MAGMA_ROOT=${OLCF_MAGMA_ROOT:="${OLCF_MAGMA_ROOT}"}
export CC=${CC:=gcc}
export FC=${FC:=gfortran}
export CXX=${CXX:=g++}
export BLAS_VENDOR=${BLAS_VENDOR:=GNU}
export BML_OPENMP=${BML_OPENMP:=yes}
export BML_OMP_OFFLOAD=${BML_OMP_OFFLOAD:=no}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export BML_TESTING=${BML_TESTING:=yes}
export BML_MAGMA=${BML_MAGMA:=yes}
export BML_COMPLEX=${BML_COMPLEX:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export EXTRA_CFLAGS=${EXTRA_CFLAGS:="-fopenmp"}
export BML_ELLBLOCK_MEMPOOL=${BML_ELLBLOCK_MEMPOOL:=no}

#use jsrun to run tests on a compute node
export BML_NONMPI_PRECOMMAND=${BML_NONMPI_PRECOMMAND:=${JSRUN}}
export BML_NONMPI_PRECOMMAND_ARGS=${BML_NONMPI_PRECOMMAND_ARGS:="-n1;-a1;-g1;-c7"}

./build.sh install
