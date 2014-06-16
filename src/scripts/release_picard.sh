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
	echo "Tags Github Picard source, checks out and builds sources, uploads build results to Sourceforge.">&2
    echo "-t <tmpdir>                Build in <tmpdir>.  Default: $TMPDIR." >&2
    echo "-u <sourceforge-user> Sourceforge username.  Default: $USERNAME." >&2
}

function tag_exists() {
    git tag | grep -q "$1$"
    if test $? = 0
        then return 0
        else return 1
    fi
}

function remote_does_not_exist() {
    git ls-remote $1 2>/dev/null 1>/dev/null
    if test $? = 0
        then return 1
        else return 0
    fi
}

function remote_tag_does_not_exist() {
    git ls-remote --tags $2 | grep -q "$1$";
    if test $? = 0
        then return 0
        else return 1
    fi
}


# This method called once for picard and once for htsjdk
function tag_it() {
    # tag must not exist
    if tag_exists $RELEASE_ID
    then echo "ERROR: Tag $RELEASE_ID locally already exists"
         return 1
    fi

    # remote must exist
    if remote_does_not_exist $REMOTE
    then echo "ERROR: Remote $REMOTE does not exist"
         return 1
    fi

    # tag at remote must not exist
    if remote_tag_does_not_exist $RELEASE_ID $REMOTE
    then echo "ERROR: Tag $RELEASE_ID at remote $REMOTE already exists"
         return 1
    fi

    # tag the branch locally then push to remote
    echo Tagging master as $RELEASE_ID and pushing the tag to $REMOTE
    # NB: we could use annotated tags in the future to store release notes, etc.
    git tag $RELEASE_ID
    git push $REMOTE $RELEASE_ID # TODO: should we check this return value in case someone made a tag since we last checked?
}

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

# Require actual Java 1.6.  This is not necessary for compiling, because can run 1.7 with -target 1.6,
# but this is necessary in order to force unit tests to run with 1.6.
(echo $JAVA_HOME | fgrep -q 1.6 ) || { echo "JAVA_HOME $JAVA_HOME is not 1.6" ; exit 1; }
java_version=`java -version 2>&1 | fgrep -i version`
(echo $java_version | fgrep -q 1.6. ) || { echo "java -version: $java_version is not 1.6"; exit 1; }

PICARDGITROOT=git@github.com:broadinstitute/picard.git
REMOTE=origin

RELEASE_ID=$1

# Since releases are lexically sorted, need to filter in order to have 1.1xx be at the bottom.
PICARD_PREV_RELEASE_ID=`git ls-remote --tags | grep -v "{}$" | awk '{print $2}' | sed -e "s_.*/__g" | egrep '[.]\d\d\d' | tail -1`

if [[ -e $TMPDIR/picard ]]
then echo "$TMPDIR/picard already exists.  Please remove or specify a different TMPDIR." >&2
        exit 1
fi
cd $TMPDIR

# clone
git clone $PICARDGITROOT picard
cd picard
ant clone-htsjdk
ant clean # clean shouldn't be necessary, but no harm

# Since releases are lexically sorted, need to filter in order to have 1.1xx be at the bottom.
PICARD_PREV_RELEASE_ID=`git ls-remote --tags | grep -v "{}$" | awk '{print $2}' | sed -e "s_.*/__g" | egrep '[.]\d\d\d' | tail -1`

HTSJDK_PREV_RELEASE_ID=$(cd htsjdk; git ls-remote --tags | grep -v "{}$" | awk '{print $2}' | sed -e "s_.*/__g" | egrep '[.]\d\d\d' | tail -1)



# Tag in both repos
for sandbox in . htsjdk
do pushd $sandbox
   tag_it || exit 1
   popd
done

ant -lib lib/ant test-htsjdk test

ant -lib lib/ant clean all javadoc

mkdir -p deploy/picard-tools/$RELEASE_ID

git log --name-status ${PICARD_PREV_RELEASE_ID}..${RELEASE_ID} > deploy/picard-tools/$RELEASE_ID/README.txt

(cd htsjdk; git log --name-status  ${HTSJDK_PREV_RELEASE_ID}..${RELEASE_ID}) >> deploy/picard-tools/$RELEASE_ID/README.txt

echo 'Edit release notes and exit editor when finished.'

$EDITOR deploy/picard-tools/$RELEASE_ID/README.txt

cp dist/picard-tools-$RELEASE_ID.zip deploy/picard-tools/$RELEASE_ID/

# Make all files to be pushed to Sourceforge writable by group so that another Picard admin can overwrite them.

chmod -R gu+rw javadoc deploy dist

find javadoc deploy dist -type d -exec chmod g+s '{}' ';' 

scp -p -r javadoc $USERNAME,picard@web.sourceforge.net:htdocs

cd deploy
scp -p -r picard-tools/$RELEASE_ID $USERNAME,picard@web.sourceforge.net:/home/frs/project/p/pi/picard/picard-tools/

cd ../dist/html
scp -p inc/*.shtml program_usage/*.shtml $USERNAME,picard@web.sourceforge.net:htdocs/inc
scp -p *.css *.shtml *.html $USERNAME,picard@web.sourceforge.net:htdocs

