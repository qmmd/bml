#!/bin/bash

set -u -e -x

echo "Starting"

env

ls -lh

bundle install
bundle exec danger

BML_OPENMP=no VERBOSE_MAKEFILE=yes EMACS=emacs26 ./build.sh check_indent
