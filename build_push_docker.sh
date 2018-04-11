#!/usr/bin/env bash

# This script is used to build and deploy docker images for Picard

if [[ "$1" == "" ]]
then
    echo "Usage: build_push_docker.sh <git-tag>"
    exit 1
fi

declare -r TAG=${1}

echo "Will build and push the following docker images:"
echo "broadinstitute/picard:${TAG}"
echo "broadinstitute/picard-cloud:${TAG}"

read -p "Is this really what you want to do? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    docker build -t broadinstitute/picard:${TAG} --build-arg build_command=shadowJar --build-arg jar_name=picard.jar .
    docker build -t us.gcr.io/broad-gotc-prod/picard-cloud:${TAG} --build-arg build_command=cloudJar --build-arg jar_name=picardcloud.jar .

    docker push broadinstitute/picard:${TAG}
    docker push us.gcr.io/broad-gotc-prod/picard-cloud:${TAG}
fi
