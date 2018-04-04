rm -r build install

# No MPI, yes GPU
CC=gcc FC=gfortran \
    BML_GPU=yes GPU_ARCH=sm_60 \
    CUDA_TOOLKIT_ROOT_DIR=/projects/opt/centos7/cuda/9.0 \
    BLAS_VENDOR=Intel CMAKE_BUILD_TYPE=Release \
    BML_OPENMP=yes \
    CMAKE_INSTALL_PREFIX=/home/smm/bml/bml/install \
    ./build.sh configure

# No MPI
#CC=gcc FC=gfortran BLAS_VENDOR=Intel CMAKE_BUILD_TYPE=Release BML_OPENMP=yes CMAKE_INSTALL_PREFIX=/home/smm/bml/bml/install ./build.sh configure

# With MPI Release
#CC=mpicc FC=mpif90 BLAS_VENDOR=Intel CMAKE_BUILD_TYPE=Release BML_OPENMP=yes BML_MPI=yes CMAKE_INSTALL_PREFIX=/home/smm/bml/bml/install ./build.sh configure

# With MPI Debug
#CC=mpicc FC=mpif90 BLAS_VENDOR=Intel CMAKE_BUILD_TYPE=Debug BML_OPENMP=yes BML_MPI=yes CMAKE_INSTALL_PREFIX=/home/smm/bml/bml/install ./build.sh configure

# With MPI Release setenv
#setenv CC mpicc;setenv FC mpif90;setenv BLAS_VENDOR Intel;setenv CMAKE_BUILD_Type Release;setenv BML_OPENMP yes;setenv BML_MPI yes;setenv CMAKE_INSTALL_PREFIX /usr/projects/infmodels/smm/qmd/bml/bml/install; setenv BML_TESTING yes; ./build.sh configure
