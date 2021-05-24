#!/bin/bash

# Make sure all the paths are correct

rm -r build
rm -r install

MY_PATH=$(pwd)

export CC=${CC:=mpicc}
export FC=${FC:=mpif90}
export CXX=${CXX:=mpic++}
export BLAS_VENDOR=${BLAS_VENDOR:=GNU}
#export BML_MPI=${BML_MPI:=yes}
export BML_MPI_NONDIST=${BML_MPI_NONDIST:=yes}
export BML_OPENMP=${BML_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export BML_TESTING=${BML_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export EXTRA_CFLAGS=${EXTRA_CFLAGS:=""}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:=""}
export MAGMA_ROOT=${MAGMA_ROOT:="${HOME}/ecp/magma-2.5.2-lib"}
#export MAGMA_ROOT=${MAGMA_ROOT:="${HOME}/ecp/magma-2.5.2"}
export BML_MAGMA=${BML_MAGMA:=yes}

./build.sh configure

                                                                                                                                                                                              
                                    
