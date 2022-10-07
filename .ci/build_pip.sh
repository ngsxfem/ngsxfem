#! /bin/bash
set -e


/opt/python/cp39-cp39/bin/python .ci/fix_auditwheel_policy.py

# 38 39

export ORIGINAL_PATH="$PATH"
for pyversion in 38 39 310
do
    #export pyversion=310
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    export NGSOLVE_VERSION=`python external_dependencies/ngsolve/tests/get_python_version_string_from_git.py external_dependencies/ngsolve`
    export PATH="$ORIGINAL_PATH:$PYDIR"
    pip install ngsolve
    #RUN $PYDIR/pip install -vvv .
    pip wheel -vvv .
    rm -rf _skbuild
    auditwheel repair xfem*.whl
    rm -f *.whl

    # avx2 build:
    pip uninstall -y ngsolve
    pip uninstall -y netgen-mesher
    
    pip install ngsolve-avx2>=$NGSOLVE_VERSION
    NETGEN_ARCH=avx2 pip wheel -vvv .
    
    auditwheel repair xfem*.whl
    rm -f *.whl
    
done


pip install -U twine
#twine upload wheelhouse/*manylinux*.whl
