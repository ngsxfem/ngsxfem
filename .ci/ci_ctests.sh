#!/bin/bash

# this script is called from the repo/basedir

set -e
echo "pwd: ${PWD}"
ls -al .

cd install
# this adjusts for different python versions and is independent of installation to site- or dist-packages
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"$(pwd)/`find -name xfem`/.."
export PYTHONPATH="$(pwd)/`find -name xfem`/..":$PYTHONPATH
cd -

echo $PYTHONPATH

cd build

if [ $1 == "demos" ]; then
  ctest -V -R 'py_demo'
fi

if [ $1 == "mayfail" ]; then
    ctest -V -R 'pymayfailtests'
fi

if [ $1 == "pytests" ]; then
    ctest -V -R 'pytests' --output-junit ctest-results.xml
fi
