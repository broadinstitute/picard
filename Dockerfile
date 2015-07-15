FROM broadinstitute/java-baseimage
MAINTAINER Broad Institute DSDE <dsde-engineering@broadinstitute.org>

COPY dist/picard.jar /usr/picard/picard.jar
WORKDIR /usr/picard

ENTRYPOINT ["java", "-jar", "picard.jar"]
