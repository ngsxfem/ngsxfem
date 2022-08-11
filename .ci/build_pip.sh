#! /bin/bash
set -e

#38 39
for pyversion in 310
do
    export pyversion=310
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    export PATH="$PATH:$PYDIR"
    $PYDIR/pip install ngsolve
    #RUN $PYDIR/pip install -vvv .
    $PYDIR/pip wheel -vvv .
    rm -rf _skbuild
    
done
$PYDIR/python .ci/fix_auditwheel_policy.py
auditwheel repair *.whl

$PYDIR/pip install -U twine
#$PYDIR/twine upload wheelhouse/*manylinux*.whl
