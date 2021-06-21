#!/bin/bash
# Make sure all the paths are correct

#module swap PrgEnv-cray/8.0.0 PrgEnv-gnu/8.0.0
#module load cmake

rm -r build
rm -r install

MY_PATH=$(pwd)

export CC=${CC:=cc}
export FC=${FC:=ftn}
export CXX=${CXX:=CC}

#"-L$ENV{ESSL_DIR}/lib64"
#setenv("CRAY_LIBSCI_DIR","/opt/cray/pe/libsci/21.04.1.1")

export BLAS_VENDOR=${BLAS_VENDOR:=CRAY}
export BML_OPENMP=${BML_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export BML_TESTING=${BML_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export EXTRA_CFLAGS=${EXTRA_CFLAGS:=""}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:=""}

./build.sh configure

                                                                                                                                                                                              
                                    
