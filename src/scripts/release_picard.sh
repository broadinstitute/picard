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

function usage () {
    echo "USAGE: $PROGNAME <release-id>" >&2
	echo "Tags Github Picard source, checks out and builds sources, uploads build results to Sourceforge.">&2
    echo "-t <tmpdir>                Build in <tmpdir>.  Default: $TMPDIR." >&2
}

function create_release () {
	local token="$1";
	local owner="$2";
	local repo="$3";
	local tag_name="$4";
	local target_commitish="$5";
	local name="$6";
	local body="$7";
	local draft="$8";
	local prerelease="$9";

	local payload="\"tag_name\":\"$tag_name\"";
	payload="$payload,\"target_commitish\":\"$target_commitish\"";
	payload="$payload,\"name\":\"$name\"";
	payload="$payload,\"body\":\"$body\"";
	payload="$payload,\"draft\":$draft";
	payload="$payload,\"prerelease\":$prerelease";
	payload="{$payload}";

	RELEASE_RESPONSE=$(curl --fail -s -S -X POST \
		https://api.github.com/repos/$owner/$repo/releases \
		-A "create-release" \
		-H "Accept: application/vnd.github.v3+json" \
		-H "Content-Type: application/json" \
		-H "Authorization: token $token" \
		-d "$payload");

	# NB: we must set the RELEASE_GITHUB_ID as the ID in the returned json response
	export RELEASE_GITHUB_ID=$(echo "$RELEASE_RESPONSE" | sed -e 's_",.*__g' -e 's_.*/__g')
}

function upload_asset () {
	local token="$1";
	local owner="$2";
	local repo="$3";
	local name="$4";
	local content_type="$5";
	local file="$6";
	local id="$7";

	curl --fail -s -S -X POST \
		https://uploads.github.com/repos/$owner/$repo/releases/$id/assets?name=$name \
		-A "upload-asset" \
		-H "Accept: application/vnd.github.v3+json" \
		-H "Content-Type: $content_type" \
		-H "Authorization: token $token" \
		--progress-bar \
		--data-binary @"$file";
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

while getopts "ht:" options; do
  case $options in
    t ) TMPDIR=$OPTARG;;
    h ) usage;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $(($OPTIND - 1))

if [ -z $GITHUB_USER_TOKEN ]
then echo "ERROR: environment variable GITHUB_USER_TOKEN must be set." >&2
	usage
	exit 1
fi

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
GHPAGES_BRANCH="gh-pages"

RELEASE_ID=$1

# Since releases are lexically sorted, need to filter in order to have 1.1xx be at the bottom.
PICARD_PREV_RELEASE_ID=`git ls-remote --tags | grep -v "{}$" | awk '{print $2}' | sed -e "s_.*/__g" | egrep '[.]\d\d\d' | tail -1`

if [[ -e $TMPDIR/picard ]]
then echo "$TMPDIR/picard already exists.  Please remove or specify a different TMPDIR." >&2
        exit 1
fi
echo "Using TMPDIR: $TMPDIR";
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

# Move the javadoc directory to a temporary location
mv javadoc tmp_javadoc

# Copy over javadoc for htsjdk since we are in the picard directory
# NB: need to move javadoc to a tmp directory since the javadoc 
# directory in the gh-pages branch may already exist.
cd htsjdk
mkdir tmp_javadoc
cp -r ../tmp_javadoc/htsjdk tmp_javadoc/.
cd ../

# Update the javadoc
for sandbox in . htsjdk
do pushd $sandbox
	if [ "." == $sandbox ]; then
		sandbox="picard";
	fi
	echo "Updating the javadoc for $sandbox"
	# Checkout the gh-pages branch
	git checkout -b $GHPAGES_BRANCH $REMOTE/$GHPAGES_BRANCH
	# Copy over from the tmp javadoc directory
	if [ ! -d javadoc ]; then
		mkdir javadoc;
	fi
	rsync -avP --delete-after tmp_javadoc/* javadoc/.
	# Remove the tmp directory as we no longer need it
	rm -r tmp_javadoc
	# Add the new javadoc files
	find javadoc/$sandbox | xargs git add
	# Commit!
	git commit -m "Updating javadoc for release: $RELEASE_ID"
	# NB: assumes the push will not fail
	git push $REMOTE $GHPAGES_BRANCH
	# Reset the repository to master 
	git checkout master
	echo "Updated the javadoc for $sandbox"
	popd
done

# Publish a release and upload assets
echo "Creating a release on github for htsjdk and picard"
create_release $GITHUB_USER_TOKEN samtools htsjdk $RELEASE_ID "" $RELEASE_ID "Release $RELEASE_ID" "false" "false";
create_release $GITHUB_USER_TOKEN broadinstitute picard $RELEASE_ID "" $RELEASE_ID "Release $RELEASE_ID" "false" "false";
echo "Github release id: $RELEASE_GITHUB_ID"
echo "Updating the release zip and README.txt to github"
upload_asset $GITHUB_USER_TOKEN broadinstitute picard picard-tools-$RELEASE_ID.zip "application/zip" deploy/picard-tools/$RELEASE_ID/picard-tools-$RELEASE_ID.zip $RELEASE_GITHUB_ID;
upload_asset $GITHUB_USER_TOKEN broadinstitute picard README.txt "application/zip" deploy/picard-tools/$RELEASE_ID/README.txt $RELEASE_GITHUB_ID;

# Update the website
echo "Updating the website"
# Assumes the gh-pages branch is already locally created
git checkout $GHPAGES_BRANCH;
cd dist/html
cp inc/*.html program_usage/*.html picard-metric-definitions.html ../../_includes/.
cd ../../
find _includes | xargs git add
git commit -m "Adding website files for $RELEASE_ID"
git push $REMOTE $GHPAGES_BRANCH

# Move back to master just in case
git checkout master

echo "Release was successful!"
