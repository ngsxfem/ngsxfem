#!/bin/bash

set -e

# cat ~/inst/netgen/bin/ngscxx
# echo "cmake-tests"
# cd build
# export NETGENDIR="${HOME}/inst/netgen/bin"
# export PATH="${NETGENDIR}:${PATH}"
# export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${HOME}/inst/netgen/lib"
# export PYTHONPATH="${PYTHONPATH}:${HOME}/inst/netgen/lib/python3/dist-packages"
# echo $PATH
# echo $(which ngscxx)

if [ $1 == "tutorial" ]; then
  ctest -V -R 'py_tutorial'
fi


if [ $1 == "mayfail" ]; then
    #nix-shell -p blas liblapack gcc gfortran gfortran.cc.lib xorg.libXmu zlib cmake mesa_noglu icu python36 suitesparse tcl-8_5 tk-8_5 python36Packages.pytest --run
    ctest -V -R 'pymayfailtests'
fi

if [ $1 == "remaining" ]; then
    #nix-shell -p blas liblapack gcc gfortran gfortran.cc.lib xorg.libXmu zlib cmake mesa_noglu icu python36 suitesparse tcl-8_5 tk-8_5 python36Packages.pytest --run
    ctest -V -R 'pytests'
fi



if [ $1 == "go4quads-tests" ]; then
   cd ../cutint/py_demos/
   #nix-shell -p blas liblapack gcc gfortran gfortran.cc.lib xorg.libXmu zlib cmake mesa_noglu icu python36 suitesparse tcl-8_5 tk-8_5 --run "
   python3 area_of_a_circle_quads.py
fi
