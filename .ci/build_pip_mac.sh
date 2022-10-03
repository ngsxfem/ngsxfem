#! /bin/bash
set -e

export PYDIR=/usr/local/Frameworks/Python.framework/Versions/$1/bin

$PYDIR/python3 --version
$PYDIR/python3 -m venv ../venv_ngs

source ../venv_ngs/bin/activate
$PYDIR/pip3 install ngsolve twine wheel

export CMAKE_OSX_ARCHITECTURES='x86_64'
$PYDIR/python3 setup.py bdist_wheel --plat-name macosx-10.15_x86_64 -j3
$PYDIR/python3 -m twine upload  --repository testpypi dist/*.whl
