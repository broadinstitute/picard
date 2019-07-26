#!/usr/bin/env bash

# This script is used to build and deploy docker images for Picard

if [[ "$1" == "" ]]
then
    echo "Usage: build_push_docker.sh <git-tag>"
    exit 1
fi

declare -r TAG=${1}

declare -r PICARD_TAG=broadinstitute/picard:${TAG}
declare -r PICARD_CLOUD_TAG=us.gcr.io/broad-gotc-prod/picard-cloud:${TAG}

echo "Will build and push the following docker images:"
echo ${PICARD_TAG}
echo ${PICARD_CLOUD_TAG}

read -p "Is this really what you want to do? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    docker build -t ${PICARD_TAG} --build-arg build_command=shadowJar --build-arg jar_name=picard.jar .
    docker build -t ${PICARD_CLOUD_TAG} --build-arg build_command=cloudJar --build-arg jar_name=picardcloud.jar .

    docker push ${PICARD_TAG}
    gcloud docker -- push ${PICARD_CLOUD_TAG}
fi