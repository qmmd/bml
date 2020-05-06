#!/bin/bash

set -x -e

: ${IMAGE_TAG:=bml-ci}
: ${IMAGE_VERSION:=2}
: ${PUSH_IMAGE:=no}

for workflow in build lint docs; do
  docker build --tag nicolasbock/${IMAGE_TAG}-${workflow}:${IMAGE_VERSION} ci-images/${workflow}
  if [[ ${PUSH_IMAGE} = yes ]]; then
    docker push nicolasbock/${IMAGE_TAG}-${workflow}:${IMAGE_VERSION}
  fi
done
