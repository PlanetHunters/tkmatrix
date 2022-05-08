#!/bin/bash

git_tag=$1
sed -i '5s/.*/version = "'${git_tag}'"/' setup.py
git add setup.py
git commit -m "Preparing release ${git_tag}"
git tag ${git_tag} -m "New release"
git push && git push --tags

