#!/usr/bin/env bash

# Example Usage: ./docker_helper.sh -j "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m" MarkDuplicates INPUT=fake.bam ...
usage() {
  cat << EOF
Usage: $0 [-j <JVM arguments>] tool_name tool_arguments

Run Picard with JVM args and program args if present. Assumes picard.jar lives in the same directory.
JVM arguments should be a quoted, space separated list.

To run picardcloud.jar, use the '-c' argument.

Example Usage:
./docker_helper.sh -j "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m" MarkDuplicates INPUT=fake.bam ...

EOF
  exit 1
}

JAR=picard.jar
while getopts "j:ch" OPTION; do
  case $OPTION in
    j) JVM_ARGS=$OPTARG;;
    c) JAR=picardcloud.jar;;
    h) usage; exit 0;;
    [?]) usage; exit 1;;
  esac
done
shift $(expr $OPTIND - 1)
TOOL_WITH_ARGS=$@

java ${JVM_ARGS} -jar /usr/picard/${JAR} ${TOOL_WITH_ARGS}
