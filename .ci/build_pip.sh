#! /bin/bash
set -e


/opt/python/cp39-cp39/bin/python .ci/fix_auditwheel_policy.py

# 38 39

export ORIGINAL_PATH="$PATH"
for pyversion in 310
do
    #export pyversion=310
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    $PYDIR/python .ci/versions_to_files.py
    export NGSOLVE_VERSION=`cat ngsolve.version`
    export PATH="$ORIGINAL_PATH:$PYDIR"
    pip install ngsolve
    #RUN $PYDIR/pip install -vvv .
    pip wheel -vvv .
    rm -rf _skbuild

    # avx2 build:
    pip uninstall -y ngsolve
    pip uninstall -y netgen-mesher
    
    pip install ngsolve-avx2>=$NGSOLVE_VERSION
    NETGEN_ARCH=avx2 pip wheel -vvv .
    
done

auditwheel repair xfem*.whl
rm *.whl

pip install -U twine
#twine upload wheelhouse/*manylinux*.whl
