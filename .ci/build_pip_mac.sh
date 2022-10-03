#! /bin/bash
set -e

export PYDIR=$Python3_ROOT_DIR/bin
$PYDIR/python3 --version
export NGSOLVE_VERSION=`$PYDIR/python3 external_dependencies/ngsolve/tests/get_python_version_string_from_git.py external_dependencies/ngsolve`

$PYDIR/pip3 install ngsolve>=NGSOLVE_VERSION

export NETGEN_Dir=$PYDIR/../lib/python$1/site-packages/netgen/cmake
export NGSolve_Dir=$PYDIR/../lib/python$1/site-packages/ngsolve/cmake
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$NGSolve_Dir:$NETGEN_Dir
export CMAKE_OSX_ARCHITECTURES='x86_64'

# $PYDIR/python3 setup.py bdist_wheel --plat-name macosx-10.15_x86_64
# $PYDIR/python3 -m twine upload --username __token__ --repository testpypi dist/*.whl
$PYDIR/pip wheel -vvv .