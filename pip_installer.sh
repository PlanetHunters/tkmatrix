#!/bin/bash

rm tests.log
rm dist* -r
set -e
tox -r > tests.log
tests_results=$(cat tests.log | grep "congratulations")
if ! [[ -z ${tests_results} ]]; then
  git_tag=$1
  sed -i '5s/.*/version = "'${git_tag}'"/' setup.py
  git add setup.py
  git commit -m "Preparing release ${git_tag}"
  git tag ${git_tag} -m "New release"
  git push
  python3 setup.py sdist bdist_wheel
  python3 -m twine upload dist/*
else
  echo "TESTS FAILED. See tests.log"
fi
rm dist* -r
