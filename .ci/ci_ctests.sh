#!/bin/bash

set -e
echo "pwd: ${PWD}"
ls -al .
# echo "cmake-tests"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${PWD}/install/lib/python3.10/dist-packages/xfem/"
export PYTHONPATH="${PYTHONPATH}:${PWD}/install/lib/python3.10/dist-packages"
# echo "${PWD}/install/lib/python3/dist-packages"
# ls -al ${PWD}/install/lib/python3/dist-packages
# export | grep PYTHONPATH

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
