#!/bin/bash

set -e -u -x

basedir=$(readlink --canonicalize $(dirname $0)/..)

[[ -f ${basedir}/scripts/ci-defaults.sh ]] && . ${basedir}/scripts/ci-defaults.sh

export CC=gcc-11
export CXX=g++-11
export FC=gfortran-11
export BUILD_SHARED_LIBS=no
export BML_OPENMP=no
export BML_INTERNAL_BLAS=no
export TESTING_EXTRA_ARGS="-R C-.*-double_real"
export BML_VALGRIND=yes

${basedir}/build.sh --debug testing
