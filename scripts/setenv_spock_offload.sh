#!/bin/bash

module load craype-accel-amd-gfx908
module load rocm/4.1.0
module load cmake
export LD_LIBRARY_PATH="$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH"

