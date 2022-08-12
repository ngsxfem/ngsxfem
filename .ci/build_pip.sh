#! /bin/bash
set -e

# 38 39

for pyversion in 310
do
    export NGSOLVE_VERSION=`$PYDIR/python external_dependencies/ngsolve/tests/get_python_version_string_from_git.py external_dependencies/ngsolve`
    #export pyversion=310
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    export PATH="$PATH:$PYDIR"
    $PYDIR/pip install ngsolve
    #RUN $PYDIR/pip install -vvv .
    $PYDIR/pip wheel -vvv .
    rm -rf _skbuild

    # avx2 build:
    $PYDIR/pip uninstall -y netgen-mesher
    $PYDIR/pip uninstall -y ngsolve
    
    $PYDIR/pip install ngsolve-avx2==$NGSOLVE_VERSION
    NETGEN_ARCH=avx2 $PYDIR/pip wheel -vvv .
    
done
$PYDIR/python .ci/fix_auditwheel_policy.py
auditwheel repair *.whl
rm *.whl

$PYDIR/pip install -U twine
$PYDIR/twine upload wheelhouse/*manylinux*.whl
