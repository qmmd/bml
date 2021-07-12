#!/bin/bash

set -e -u -x

basedir=$(readlink --canonicalize $(dirname $0)/..)

[[ -f ${basedir}/scripts/ci-defaults.sh ]] && . ${basedir}/scripts/ci-defaults.sh

export CC=gcc-5
export CXX=g++-5
export FC=gfortran-5
export BUILD_SHARED_LIBS=yes
export BML_OPENMP=yes
export BML_INTERNAL_BLAS=yes

${basedir}/build.sh testing
