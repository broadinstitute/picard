ARG BASE_DOCKER=eclipse-temurin:17-jdk
ARG release=false

FROM ${BASE_DOCKER} as build
LABEL stage=buildStage

# Install git for building
RUN apt-get update && \
    apt-get --no-install-recommends install -y \
        git

# Assumes Dockerfile lives in root of the git repo. Pull source files into container
COPY / /usr/picard/
WORKDIR /usr/picard

# download gradle then build
RUN ./gradlew -Drelease=${release} \
   clean \
   printVersion \
   shadowJar

FROM ${BASE_DOCKER} as final
MAINTAINER Broad Institute DSDE <dsde-engineering@broadinstitute.org>

# Install R
RUN apt-get update && \
    apt-get --no-install-recommends install -y \
        r-base &&\
    apt-get clean autoclean && \
    apt-get autoremove -y

RUN mkdir /usr/picard/

COPY --from=build /usr/picard/build/libs/picard.jar /usr/picard/

RUN mkdir /usr/working
WORKDIR /usr/working
