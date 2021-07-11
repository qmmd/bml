#!/bin/bash

set -e -u -x

basedir=$(readlink --canonicalize $(dirname $0)/..)

[[ -f ${basedir}/scripts/ci-defaults.sh ]] && . ${basedir}/scripts/ci-defaults.sh

export CC=gcc-4.8
export CXX=g++-4.8
export FC=gfortran-4.8
export BUILD_SHARED_LIBS=yes
export BML_OPENMP=no
export BML_INTERNAL_BLAS=no

${basedir}/build.sh testing
