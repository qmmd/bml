#!/bin/bash

set -e -u -x

basedir=$(readlink --canonicalize $(dirname $0)/..)

[[ -f ${basedir}/scripts/ci-defaults.sh ]] \
    && . ${basedir}/scripts/ci-defaults.sh

export CC=mpicc
export CXX=mpic++
export FC=mpifort
export BUILD_SHARED_LIBS=no
export BML_OPENMP=no
export BML_INTERNAL_BLAS=no
export BML_MPI=yes
export TESTING_EXTRA_ARGS="-R MPI-C-.*-single_real"
export EXTRA_LINK_FLAGS=-lscalapack-openmpi
export BML_SCALAPACK=yes

${basedir}/build.sh --debug testing
