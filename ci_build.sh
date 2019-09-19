#!/bin/bash
mkdir -p build
mkdir -p ~/ngsxfem-reports
mkdir -p ~/inst
mkdir -p ~/inst/netgen
touch ~/ngsxfem-reports/lastngs.version
chmod 755 ~/ngsxfem-reports/lastngs.version
git submodule update --init --recursive
export NGS_VERSION=`git submodule | awk '{print $1;}'`
export NGS_INSTALLED_VERSION=`cat ~/ngsxfem-reports/lastngs.version | awk '{print $1;}'`
echo NGS_VERSION is $NGS_VERSION
echo NGS_INSTALLED_VERSION is $NGS_INSTALLED_VERSION
export BUILD_NGS="ON"
test "$NGS_VERSION" = "$NGS_INSTALLED_VERSION" && echo "versions match" && [ -d ~/ngsxfem-reports/last-inst ] && export BUILD_NGS="OFF" && cp -rv ~/ngsxfem-reports/last-inst/* ~/inst/netgen/
test "$NGS_VERSION" != "$NGS_INSTALLED_VERSION" && echo "versions do not match"
cd build
cmake -DCMAKE_INSTALL_PREFIX=~/inst/netgen -DCMAKE_BUILD_TYPE=RELEASE -DUSE_CCACHE=ON -DUSE_GUI=OFF -DTCL_INCLUDE_PATH=/usr/include/tcl8.5/ -DBUILD_NGSOLVE=$BUILD_NGS -DBUILD_NGSOLVE_THREADS=15 -GNinja ../.
ninja install
echo $NGS_VERSION > ~/ngsxfem-reports/lastngs.version
rm -rf ~/ngsxfem-reports/last-inst
cp -r ~/inst/netgen ~/ngsxfem-reports/last-inst
chmod 755 ~/ngsxfem-reports/lastngs.version
