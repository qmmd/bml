#!/bin/bash

module restore PrgEnv-cray
module unload cce
module use /home/groups/coegroup/share/coe/modulefiles
module load hipmagma
module load cmake
export LD_LIBRARY_PATH="$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH"

