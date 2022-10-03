#! /bin/bash
set -e

$Python3_ROOT_DIR/python3 --version
$Python3_ROOT_DIR/python3 -m venv ../venv_ngs

source ../venv_ngs/bin/activate
$Python3_ROOT_DIR/pip3 install ngsolve twine wheel

export CMAKE_OSX_ARCHITECTURES='x86_64'
$Python3_ROOT_DIR/python3 setup.py bdist_wheel --plat-name macosx-10.15_x86_64 -j3
$Python3_ROOT_DIR/python3 -m twine upload  --repository testpypi dist/*.whl
