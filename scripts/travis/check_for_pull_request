#!/usr/bin/env bash

#query github for information about a specific pull request
curl -v -H "Authorization: token $GITHUB_API_TOKEN" https://api.github.com/repos/broadinstitute/picard/pulls\?state\=open\&head="broadinstitute:${TRAVIS_BRANCH}"
