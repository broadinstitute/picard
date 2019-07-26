FROM marketplace.gcr.io/google/debian9
MAINTAINER Broad Institute DSDE <dsde-engineering@broadinstitute.org>

ARG build_command=shadowJar
ARG jar_name=picard.jar

# Install ant, git for building
RUN apt-get update && \
    apt-get --no-install-recommends install -y --force-yes \
        git \
        r-base \
        ant \
        openjdk-8-jdk && \
    apt-get clean autoclean && \
    apt-get autoremove -y

# Assumes Dockerfile lives in root of the git repo. Pull source files into container
COPY / /usr/picard/
WORKDIR /usr/picard

# Build the distribution jar, clean up everything else
RUN ./gradlew ${build_command} && \
    mv build/libs/${jar_name} picard.jar && \
    ./gradlew clean && \
    rm -rf src && \
    rm -rf gradle && \
    rm -rf .git && \
    rm gradlew && \
    rm build.gradle

RUN mkdir /usr/working
WORKDIR /usr/working