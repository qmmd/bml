#!/bin/bash

module load craype-accel-amd-gfx908
module unload cray-mvapich2 cray-libsci
module use /home/groups/coegroup/share/coe/modulefiles
module load rocm/4.1.1
module swap cce cce/11.0.4
module load ompi/4.1.0/cce/11.0.4/rocm/4.1.1
module load blas
#module load rocm/4.2.0
module load cmake
export LD_LIBRARY_PATH="$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH"

