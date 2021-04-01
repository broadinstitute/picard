#!/usr/bin/env bash

#query github for information about a specific pull request
PULL_REQUESTS=$(curl -v -H "Authorization: token $GITHUB_API_TOKEN" https://api.github.com/repos/broadinstitute/picard/pulls\?state\=open\&head="broadinstitute:${TRAVIS_BRANCH}")
if [[ $( grep -c "commits" <<< ${PULL_REQUESTS} ) -gt 0 ]]; then
  echo "pull requests non-empty"
  return 0
else
  echo "pull requests empty"
  return 1
fi