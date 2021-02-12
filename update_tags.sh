#!/bin/bash

set -x

basedir=$(readlink --canonicalize $(dirname $0))
pushd "${basedir}" || exit
readarray -t files < <(git ls-files '*.c' '*.h' '*.F90')
ctags "${files[@]}"
etags "${files[@]}"
popd || exit
