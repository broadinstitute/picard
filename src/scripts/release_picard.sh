#! /bin/bash

# The MIT License
#
# Copyright (c) $today.year The Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN

PROGNAME=`basename $0`
USERNAME=alecw

function usage () {
    echo "USAGE: $PROGNAME <release-id>" >&2
    echo "Branches Sourceforge Picard source, checks out and builds sources, uploads build results to Sourceforge.">&2
    echo "-t <tmpdir>                Build in <tmpdir>.  Default: $TMPDIR." >&2
    echo "-u <sourceforge-user> Sourceforge username.  Default: $USERNAME." >&2
}

function branch_exists() {
    if svn info $1 2>&1 | fgrep -q 'Not a valid URL'
        then return 1
        else return 0
    fi
}

SVNROOT=https://picard.svn.sourceforge.net/svnroot/picard

set -e

while getopts "ht:u:" options; do
  case $options in
    u ) USERNAME=$OPTARG;;
    t ) TMPDIR=$OPTARG;;
    h ) usage;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;

  esac
done
shift $(($OPTIND - 1))

if (( $# != 1 ))
 then echo "ERROR: Incorrect number of arguments." >&2
      usage
      exit 1
fi

if [[ x"$EDITOR" == x ]]
then echo "EDITOR environment variable must be set." >&2
       exit 1
fi

RELEASE_ID=$1

PREV_RELEASE_ID=`svn ls $SVNROOT/tags | tail -1 | sed 's/\/$//'`

if branch_exists $SVNROOT/branches/$RELEASE_ID
then echo "ERROR: $SVNROOT/branches/$RELEASE_ID already exists.">&2
       exit 1
fi

if branch_exists $SVNROOT/tags/$RELEASE_ID
then echo "ERROR: $SVNROOT/tags/$RELEASE_ID already exists.">&2
       exit 1
fi

if [[ -e $TMPDIR/Picard-public ]]
then echo "$TMPDIR/Picard-public already exists.  Please remove or specify a different TMPDIR." >&2
        exit 1
fi

svn copy -m "Release $RELEASE_ID" $SVNROOT/trunk $SVNROOT/branches/$RELEASE_ID
svn copy -m "Release $RELEASE_ID" $SVNROOT/trunk $SVNROOT/tags/$RELEASE_ID

cd $TMPDIR

mkdir Picard-public
cd Picard-public
svn co $SVNROOT/tags/$RELEASE_ID .

ant test

ant clean all javadoc

REVISION=`svn info $SVNROOT/tags/$RELEASE_ID | egrep '^Last Changed Rev: ' | awk '{print $4}'`
PREV_REVISION=`svn info $SVNROOT/tags/$PREV_RELEASE_ID | egrep '^Last Changed Rev: ' | awk '{print $4}'`

mkdir -p deploy/picard-tools/$RELEASE_ID

svn log -r $PREV_REVISION:$REVISION -v > deploy/picard-tools/$RELEASE_ID/release_notes.txt

echo 'Edit release notes and exit editor when finished.'

$EDITOR deploy/picard-tools/$RELEASE_ID/release_notes.txt

cp dist/picard-tools-$RELEASE_ID.zip deploy/picard-tools/$RELEASE_ID/
mkdir -p deploy/sam-jdk/$RELEASE_ID
cp dist/sam-$RELEASE_ID.jar deploy/sam-jdk/$RELEASE_ID/

# Make all files to be pushed to Sourceforge writable by group so that another Picard admin can overwrite them.

chmod -R gu+rw javadoc deploy dist

scp -p -r javadoc $USERNAME,picard@web.sourceforge.net:htdocs

cd deploy
scp -p -r picard-tools/$RELEASE_ID $USERNAME,picard@web.sourceforge.net:/home/frs/project/p/pi/picard/picard-tools/
scp -p -r sam-jdk/$RELEASE_ID $USERNAME,picard@web.sourceforge.net:/home/frs/project/p/pi/picard/sam-jdk/

cd ../dist/html
scp -p *.shtml program_usage/*.shtml $USERNAME,picard@web.sourceforge.net:htdocs/inc

