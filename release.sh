#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
rm tests.log
rm dist* -r
rm -r .tox
rm -r .pytest_cache
rm -r build
rm -r tkmatrix-reqs
rm -R *egg-info
conda remove -n tkmatrix-reqs --all -y
set -e
git_tag=$1
echo "GIT TAG IS " ${git_tag}
tox -r -e py310 > tests.log
tests_results=$(cat tests.log | grep "congratulations")
if ! [[ -z ${tests_results} ]]; then
  set +e
  rm tests.log
  rm -r .tox
  rm -r .pytest_cache
  rm -r build
  rm -r tkmatrix-reqs
  rm -R *egg-info
  conda remove -n tkmatrix-reqs --all -y
  set -e
  conda create -n tkmatrix-reqs python=3.10 -y
  conda activate tkmatrix-reqs
  python3 -m pip install pip -U
  python3 -m pip install numpy==1.23.5
  sed -i '5s/.*/version = "'${git_tag}'"/' setup.py
  python3 -m pip install -e .
  python3 -m pip list --format=freeze > requirements.txt
  conda deactivate
  git add requirements.txt
  git add setup.py
  git commit -m "Preparing release ${git_tag}"
  git tag ${git_tag} -m "New release"
  git push && git push --tags
else
  echo "TESTS FAILED. See tests.log"
fi
set +e
rm -R tkmatrix-reqs
rm dist* -r
rm -r .tox
rm -r .pytest_cache
rm -r build
rm -R *egg-info
conda remove -n tkmatrix-reqs --all -y
