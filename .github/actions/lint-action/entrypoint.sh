#!/bin/bash

set -u -e -x

echo "Starting"

env

ls -lh

BML_OPENMP=no VERBOSE_MAKEFILE=yes EMACS=emacs26 ./build.sh check_indent
