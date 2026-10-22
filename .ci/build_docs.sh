#!/bin/bash
# Builds the sphinx documentation into ./public
# (used by the GitLab job "pages" and the GitHub workflow "pages").
#
# Needs
#  - an installed xfem in ./install (e.g. from the job build-on-pip-ngs),
#  - ngsolve + the python packages of doc/sphinx/requirements.txt (+ ipykernel,
#    webgui_jupyter_widgets, scipy, matplotlib for the notebooks),
#  - the submodule external_dependencies/jupyter.
# Optional: literature.pdf / feature-details.pdf in the repository root
# (job paper-appendix); without them the pdf link on the literature page is removed.
set -e
export TZ=UTC
cd "$(dirname "$0")/.."

XFEM_INIT=$(find install/lib -path '*/xfem/__init__.py' | head -n 1)
[ -n "$XFEM_INIT" ] || { echo "no installed xfem found in install/lib" >&2; exit 1; }
XFEM_PATH=$(realpath "$(dirname "$XFEM_INIT")/..")
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${XFEM_PATH}
export PYTHONPATH=${XFEM_PATH}:${PYTHONPATH}
python3 -c "import ngsolve, xfem; print('ngsolve', ngsolve.__version__, '| xfem from', xfem.__file__)"

rm -rf html public

# generate html-files and copy some stuff
sphinx-apidoc -e -f -o doc/sphinx/xfem_doc $XFEM_PATH/xfem
cp external_dependencies/jupyter/{aggregation,basics,intlset,cutfem,spacetime_fitted,spacetime_unfitted,mlset_basic,mlset_pde,unfmixed}.ipynb doc/sphinx/jupyter
cp external_dependencies/jupyter/helper.py doc/sphinx/jupyter
cp -r external_dependencies/jupyter/graphics doc/sphinx/jupyter
cp demos/*.py doc/sphinx/xfem_misc/
mkdir -p doc/sphinx/xfem_misc/aggregates doc/sphinx/xfem_misc/mpi doc/sphinx/xfem_misc/spacetime doc/sphinx/xfem_misc/spdes
cp demos/aggregates/*.py doc/sphinx/xfem_misc/aggregates
cp demos/mpi/*.py doc/sphinx/xfem_misc/mpi
cp demos/spacetime/*.py doc/sphinx/xfem_misc/spacetime
cp demos/spdes/*.py doc/sphinx/xfem_misc/spdes
cp -r doc/graphics doc/sphinx/xfem_misc/
sphinx-build doc/sphinx html
cp doc/sphinx/bullets.css html/_static

# copy missing files
mkdir -p html/doc html/xfem_misc/python
cp doc/feature-details.md html/doc/feature-details.md
# compiled pdfs (job paper-appendix), available as
#   <pages-url>/doc/literature.pdf and <pages-url>/doc/feature-details.pdf
for f in literature feature-details; do
  cp $f.pdf html/doc/ 2>/dev/null || echo "no compiled $f.pdf available"
done
# no dead link on the literature page if its pdf is missing
[ -f html/doc/literature.pdf ] || sed -i '/class="pub-pdf"/d' html/xfem_misc/literature.html
cp INSTALLATION.md html/INSTALLATION.md
cp python/ngs_check.py html/xfem_misc/python/ngs_check.py

rm -f html/.doctrees/nbsphinx/jupyter/*.ipynb
rm -f html/jupyter/*.ipynb
mv html public
