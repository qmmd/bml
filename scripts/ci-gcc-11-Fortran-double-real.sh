#!/bin/bash

set -e -u -x

basedir=$(readlink --canonicalize $(dirname $0)/..)

[[ -f ${basedir}/scripts/ci-defaults.sh ]] \
    && . ${basedir}/scripts/ci-defaults.sh

export CC=gcc-11
export CXX=g++-11
export FC=gfortran-11
export BUILD_SHARED_LIBS=no
export BML_OPENMP=no
export BML_INTERNAL_BLAS=no
export TESTING_EXTRA_ARGS="-R fortran-.*-double_real"
export BML_VALGRIND=yes

# [nicolasbock] Valgrind cannot handle a Debug build with valgrind and
# fails with:
#
# parse DIE(readdwarf3.c:3619): confused by:
#  <2><60e>: Abbrev Number: 25 (DW_TAG_subrange_type)
#      DW_AT_type        : <42>
#      DW_AT_upper_bound : 4 byte block: 91 c0 7b 6
# parse_type_DIE:
# --18653-- WARNING: Serious error when reading debug info
# --18653-- When reading debug info from /home/runner/work/bml/bml/build/tests/Fortran-tests/bml-testf:
# --18653-- confused by the above DIE
#
# Bionic ships with valgrind-3.13.0. This issue might be fixed with
# later versions.
#
# Re-review once we have a backport of valgrind for Bionic.
export CMAKE_BUILD_TYPE=Release

${basedir}/build.sh --debug testing
