#!/bin/bash

basedir=$(dirname $0)
readarray -t files < <(git ls-files '*.c' '*.h' '*.F90')
pushd "${basedir}" || exit
ctags "${files[@]}"
etags "${files[@]}"
popd || exit
