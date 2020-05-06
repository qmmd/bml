#!/bin/bash

set -x

: ${IMAGE_TAG:=bml-ci}
: ${IMAGE_VERSION:=1}

for workflow in build lint docs; do
  docker build --tag nicolasbock/${IMAGE_TAG}-${workflow}:${IMAGE_VERSION} ci-images/${workflow}
  docker push nicolasbock/${IMAGE_TAG}-${workflow}:${IMAGE_VERSION}
done
