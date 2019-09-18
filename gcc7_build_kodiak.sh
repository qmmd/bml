#!/bin/bash

# Make sure all the paths are correct

module purge
module load gcc/7.4.0 cmake cudatoolkit/8.0 mkl/11.4.1

rm -r build
rm -r install

MY_PATH=$(pwd)

export CC=${CC:=gcc}
export FC=${FC:=gfortran}
export CXX=${CXX:=g++}
export BLAS_VENDOR=${BLAS_VENDOR:=MKL}
export BML_OPENMP=${BML_OPENMP:=yes}
export BML_GPU=${BML_GPU:=no}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export BML_TESTING=${BML_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export EXTRA_CFLAGS=${EXTRA_CFLAGS:=""}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:=""}
export CMAKE_Fortran_FLAGS="-fopenmp"
export CMAKE_C_FLAGS="-O3 -fopenmp"
export CUDA_TOOLKIT_ROOT_DIR=${CUDATOOLKIT_ROOT}
./build.sh install

                                                                                                                                                                                              
                                    
