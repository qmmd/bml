#!/bin/bash

set -ex

for i in $(seq 5); do
    echo "Attempt ${i} to install packages"
    if ! sudo apt-get update; then
        continue
    fi
    if sudo apt-get install --no-install-recommends ${DEVEL_PACKAGES} ${CC} ${CXX} ${FC} ${packages}; then
        exit 0
    fi
done

exit 1
