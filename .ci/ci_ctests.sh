#!/bin/bash

set -e
echo "pwd: ${PWD}"
ls -al .
echo "cmake-tests"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${PWD}/install/lib/python3/dist-packages/xfem/"
export PYTHONPATH="${PYTHONPATH}:${PWD}/install/lib/python3/dist-packages"
echo "${PWD}/install/lib/python3/dist-packages"
ls -al ${PWD}/install/lib/python3/dist-packages
export | grep PYTHONPATH

cd build

if [ $1 == "tutorial" ]; then
  ctest -V -R 'py_tutorial'
fi

if [ $1 == "mayfail" ]; then
    ctest -V -R 'pymayfailtests'
fi

if [ $1 == "remaining" ]; then
    ctest -V -R 'pytests'
fi

if [ $1 == "go4quads-tests" ]; then
   cd ../cutint/py_demos/
   python3 area_of_a_circle_quads.py
fi
