#CC=pgcc FC=pgf90 BLAS_VENDOR=Custom BML_INTERNAL_BLAS=no CUSTOM_BLAS_DIR=/home/mewall/pgi/lapack-release CUSTOM_CLAPACK_DIR=/home/mewall/pgi/clapack/install/lib CMAKE_BUILD_TYPE=Release BML_OPENMP=yes CMAKE_INSTALL_PREFIX=/home/mewall/power/bml/install ./build.sh install
CC=pgcc FC=pgf90 BLAS_VENDOR=Custom BML_INTERNAL_BLAS=no CUSTOM_BLAS_DIR=/home/mewall/pgi/lapack-release CUSTOM_CLAPACK_DIR=/home/mewall/pgi/clapack/install/lib CMAKE_BUILD_TYPE=Debug BML_OPENMP=yes CMAKE_INSTALL_PREFIX=/home/mewall/power/bml/install ./build.sh install

