#! /bin/bash
set -e

export PYDIR=$Python3_ROOT_DIR/bin
$PYDIR/python3 --version
$PYDIR/pip3 install ngsolve twine wheel pybind11-stubgen


export NETGEN_Dir=$PYDIR/../lib/python$1/site-packages/netgen/cmake
export NGSolve_Dir=$PYDIR/../lib/python$1/site-packages/ngsolve/cmake
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$NGSolve_Dir:$NETGEN_Dir
export CMAKE_OSX_ARCHITECTURES='x86_64'

cd $NGSolve_Dir
ls -lta

# $PYDIR/python setup.py bdist_wheel --plat-name macosx-10.15_x86_64
# $PYDIR/python -m twine upload --username __token__ --repository testpypi dist/*.whl
