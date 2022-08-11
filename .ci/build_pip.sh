#! /bin/bash
set -e

#38 39
for pyversion in 310

    export pyversion=310
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    $PYDIR/pip install ngsolve
    #RUN $PYDIR/pip install -vvv .
    $PYDIR/pip wheel -vvv .

    rm -rf _skbuild
$PYDIR/python fix_auditwheel_policy.py
auditwheel repair *.whl

$PYDIR/pip install -U twine
$PYDIR/twine upload wheelhouse/*manylinux*.whl