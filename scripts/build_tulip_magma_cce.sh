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
export CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS:="-L${CRAY_LIBSCI_PREFIX}/lib -L${HIP_PATH}/lib -L${MAGMA_HOME}/lib -L${ROCM_PATH}/lib -lsci_cray -lmagma -D__HIP_PLATFORM_HCC__ -DHAVE_HIP -ef -DCRAY_SDK"}
export EXTRA_CFLAGS=${EXTRA_CFLAGS:="-I${CRAY_LIBSCI_PREFIX}/include -I${HIP_PATH}/include -I${MAGMA_HOME}/include -I${ROCM_PATH}/include -DHAVE_HIP -D__HIP_PLATFORM_HCC__ -DCRAY_SDK"}
#export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-L${HIP_PATH}/lib -L/home/users/coe0209/packages/openblas -lamdhip64 -lopenblas -lm"}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-L${CRAY_LIBSCI_PREFIX}/lib -L${HIP_PATH}/lib -L${MAGMA_HOME}/lib -L${ROCM_PATH}/lib -lsci_cray -lmagma"}

./build.sh configure

pushd build
make -j8
make install
popd
