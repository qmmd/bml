#!/bin/bash

set -x

: ${IMAGE_TAG:=bml-ci}
: ${IMAGE_VERSION:=1.0}

docker build --tag nicolasbock/${IMAGE_TAG}:${IMAGE_VERSION} ci-images
docker push nicolasbock/${IMAGE_TAG}:${IMAGE_VERSION}
