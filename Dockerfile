FROM broadinstitute/java-baseimage
MAINTAINER Broad Institute DSDE <dsde-engineering@broadinstitute.org>

# Install ant, git for building
RUN apt-get update && \
    apt-get --no-install-recommends install -y --force-yes \
        git \
        ant && \
    apt-get clean autoclean && \
    apt-get autoremove -y

# Assumes Dockerfile lives in root of the git repo. Pull source files into container
COPY / /usr/picard/
WORKDIR /usr/picard

# Build the distribution jar, clean up everything else
RUN ./gradlew shadowJar && \
    mv build/libs/picard.jar picard.jar && \
    mv src/main/resources/picard/docker_helper.sh docker_helper.sh && \
    ./gradlew clean && \
    rm -rf src && \
    rm -rf gradle && \
    rm -rf .git && \
    rm gradlew && \
    rm build.gradle

RUN mkdir /usr/working
WORKDIR /usr/working

ENTRYPOINT ["/usr/picard/docker_helper.sh"]
CMD [""]
