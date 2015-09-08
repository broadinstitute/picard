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
COPY build.xml /usr/picard/build.xml
COPY src /usr/picard/src
COPY lib /usr/picard/lib
WORKDIR /usr/picard

# Clone out htsjdk. First turn off git ssl verification
RUN git config --global http.sslVerify false && git clone https://github.com/samtools/htsjdk.git

# Build the distribution jar, clean up everything else
RUN ant clean all && \
    mv dist/picard.jar picard.jar && \
    mv src/scripts/picard/docker_helper.sh docker_helper.sh && \
    ant clean && \
    rm -rf htsjdk && \
    rm -rf src && \
    rm -rf lib && \
    rm build.xml

RUN mkdir /usr/working
WORKDIR /usr/working

ENTRYPOINT ["/usr/picard/docker_helper.sh"]
CMD [""]
