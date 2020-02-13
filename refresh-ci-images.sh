#!/bin/bash

set -x

docker build --tag nicolasbock/bml-ci ci-images
docker push nicolasbock/bml-ci
