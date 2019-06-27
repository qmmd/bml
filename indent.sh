#!/bin/bash

set -x

INDENT_ARGS="-gnu -nut -i4 -bli0 -cli4 -ppi0 -cbi0 -npcs -bfda"

declare -a C_FILES
declare -a FORTRAN_FILES

if ! git status > /dev/null; then
  cat <<EOF
We are not inside a git repository. Please specify the files to indent
manually.
EOF
  exit 1
fi

if (( $# > 0 )); then
  for file in "$@"; do
    case "${file##*.}" in
      "c" | "h")
        C_FILES[${#C_FILES[@]}]="$file"
        ;;
      "F90")
        FORTRAN_FILES[${#FORTRAN_FILES[@]}]="$file"
        ;;
      *)
        echo "unknown suffix"
        ;;
    esac
  done
else
  readarray -t C_FILES < <(git ls-files -- *.c)
  readarray -t FORTRAN_FILES < <(git ls-files -- *.F90)
fi

for f in "${C_FILES[@]}" "${FORTRAN_FILES[@]}"; do
  sed -i -e 's:\s\+$::' "${f}"
done

for f in "${C_FILES[@]}"; do
  indent ${INDENT_ARGS} "${f}"
done

for f in "${FORTRAN_FILES[@]}"; do
  emacs26 --batch \
    "${f}" \
    --eval "(describe-variable major-mode)" \
    --eval "(indent-region (minibuffer-prompt-end) (point-max) nil)" \
    -f save-buffer
done
