#! /bin/bash
set -e
mkdir dist

export PYDIR=$Python3_ROOT_DIR/bin
$PYDIR/python3 --version
export NGSOLVE_VERSION=`$PYDIR/python3 external_dependencies/ngsolve/tests/get_python_version_string_from_git.py external_dependencies/ngsolve`
export NETGEN_VERSION=`$PYDIR/python3 external_dependencies/ngsolve/tests/get_python_version_string_from_git.py external_dependencies/ngsolve/external_dependencies/netgen`

$PYDIR/pip3 install pybind11-stubgen
$PYDIR/pip3 install netgen-mesher>=$NETGEN_VERSION --pre
$PYDIR/pip3 install ngsolve>=$NGSOLVE_VERSION --pre
$PYDIR/pip3 install scikit-build wheel

export PATH=/Applications/CMake.app/Contents/bin:$PATH
export NETGEN_Dir=$PYDIR/../lib/python$1/site-packages/netgen/cmake
export NGSolve_Dir=$PYDIR/../lib/python$1/site-packages/ngsolve/cmake
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$NGSolve_Dir:$NETGEN_Dir
export CMAKE_OSX_ARCHITECTURES='arm64;x86_64'

$PYDIR/pip wheel -vvv -w dist .

$PYDIR/python3 setup.py bdist_wheel --plat-name macosx-10.15-universal2 -d wheelhouse

#cd dist; rm -f netgen*.whl ngsolve*.whl
#wheel tags --platform-tag macosx_10_15_x86_64 xfem-*.whl
