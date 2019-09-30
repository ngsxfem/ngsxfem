#!/bin/bash

cd ~/inst/netgen/share/ngsxfem/report
mkdir -p ~/ngsxfem-reports
touch ~/ngsxfem-reports/log
echo "" >> ~/ngsxfem-reports/log
date >> ~/ngsxfem-reports/log
echo pipeline-id $CI_PIPELINE_ID build $CI_BUILD_ID $CI_BUILD_TAG $CI_BUILD_REF >> ~/ngsxfem-reports/log
export NETGENDIR="${HOME}/inst/netgen/bin"
export PATH="${NETGENDIR}:${PATH}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${HOME}/inst/netgen/lib"
export PYTHONPATH="${PYTHONPATH}:${HOME}/inst/netgen/lib/python3/dist-packages"
mkdir ~/ngsxfem-reports/$CI_BUILD_REF_NAME/ -p
#nix-shell -p blas liblapack gcc gfortran gfortran.cc.lib xorg.libXmu zlib cmake mesa_noglu icu python36 suitesparse tcl-8_5 tk-8_5 ninja --run
python3 ngsxfem_report.py ~/ngsxfem-reports/$CI_BUILD_REF_NAME/ $CI_PIPELINE_ID
