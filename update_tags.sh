#!/bin/bash


files=$(git ls-tree --full-tree -r HEAD \
    | grep '.\(c\|h\|F90\)$' \
    | awk '{print $4}')
ctags --recurse --C-kinds=+lxzLp --Fortran-kinds=+LP ${files}
etags ${files}
