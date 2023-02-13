#!/bin/bash

rm tests.log
rm dist* -r
rm -r .tox
rm -r .pytest_cache
rm -r build
rm -r tkmatrix-reqs
rm -R *egg-info
set -e

git_tag=$1
echo "GIT TAG IS " ${git_tag}
tox -r -e py38,py39 > tests.log
tests_results=$(cat tests.log | grep "congratulations")
if ! [[ -z ${tests_results} ]]; then
  set +e
  rm tests.log
  rm -r .tox
  rm -r .pytest_cache
  rm -r build
  rm -r tkmatrix-reqs
  rm -R *egg-info
  set -e
  python3.8 -m venv tkmatrix-reqs
  source tkmatrix-reqs/bin/activate
  python3.8 -m pip install pip -U
  python3.8 -m pip install numpy==1.22.4
  sed -i '5s/.*/version = "'${git_tag}'"/' setup.py
  python3.8 -m pip install -e .
  python3.8 -m pip list --format=freeze > requirements.txt
  deactivate
  git add requirements.txt
  git add setup.py
  git commit -m "Preparing release ${git_tag}"
  git tag ${git_tag} -m "New release"
  git push && git push --tags
else
  echo "TESTS FAILED. See tests.log"
fi
set +e
rm -R sherlockpipe-reqs
rm dist* -r
rm -r .tox
rm -r .pytest_cache
rm -r build
rm -R *egg-info
