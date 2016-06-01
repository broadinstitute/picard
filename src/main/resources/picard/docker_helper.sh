#!/usr/bin/env bash

# Example Usage: ./docker_helper.sh -j "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m" MarkDuplicates INPUT=fake.bam ...
usage() {
  cat << EOF
Usage: $0 [-j <JVM arguments>] tool_name tool_arguments

Run Picard with JVM args and program args if present. Assumes picard.jar lives in the same directory.
JVM arguments should be a quoted, space separated list.

Example Usage:
./docker_helper.sh -j "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m" MarkDuplicates INPUT=fake.bam ...

EOF
  exit 1
}
while getopts "j:h" OPTION; do
  case $OPTION in
    j) JVM_ARGS=$OPTARG;;
    h) usage; exit 0;;
    [?]) usage; exit 1;;
  esac
done
shift $(expr $OPTIND - 1)
TOOL_WITH_ARGS=$@

java ${JVM_ARGS} -jar /usr/picard/picard.jar ${TOOL_WITH_ARGS}
