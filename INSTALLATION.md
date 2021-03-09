# Installation of `ngsxfem`
We provide two (and a half) ways to setup `ngsxfem` on your machine:
* [Building/Installing through `pip`](#installation-through-pip-pip3)
* [Building/Installing from sources](#installation-from-source)
* [Running a docker image](#docker-container)

Below we discuss these installation steps in more detail. If you observe any problems with the installation, you can contact us through the [github issue tracker](https://github.com/ngsxfem/ngsxfem/issues) or the [`'NGSolve' user forum`](https://ngsolve.org/forum/index).

# Installation through `pip` (`pip3`)
Table of contents for `pip` install:
  * [Releases](#releases)
  * [Development version](#development-version)
  * [Troubleshooting](#troubleshooting)
---

We automatically publish a *source distribution* of the releases of `ngsxfem` on PyPI for installation via `pip3`. Since this is a source code distribution and not a binary distribution, the installation via `pip` requires a sufficiently new pre-installed `NGSolve` and the same prerequisites as installing `NGSolve` from source. For details of these, please consult the [corresponding section on the `NGSolve` homepage](https://ngsolve.org/docu/latest/install/installlinux.html). The necessary `NGSolve` version for compatibility will be checked for during installation and can be found in the first line of [`python/ngs_check.py`](python/ngs_check.py). Note that the version of the submodule `ngsolve` in `external dependencies` is also sufficient.

## Releases
Installation of `ngsxfem` releases through `pip` is carried with the command
``` {.shell}
pip3 install xfem
```
You may add standard `pip` options such as `--user`, `--upgrade` and/or `--verbose` or specify a concrete version, e.g. by replacing `xfem` with `xfem==1.4.2104`.
Note that the installation requires some compilations and may hence take a few minutes.

## Development version

You can also directly install (through `pip`) the latest version of ngsxfem which is the development version in the master branch on github using
``` {.shell}
pip3 install git+https://github.com/ngsxfem/ngsxfem.git@master
```

## Troubleshooting
The `pip`-installation builds `ngsxfem` from source using default parameters for the build. If you meet problems with the `pip`-installation you may want to try an [installation from source](#building-from-source)
as it gives you more fine-grained control on the installation. 

# Installation from source
You can build `ngsxfem` also directly from sources. This should work for both the releases and the master branch. As `ngsxfem` is an Add-on to `NGSolve` you require a proper version of `NGSolve`, the installation of which we discuss first before we explain the further build steps.

Table of contents for install from source:
* [Prerequisite `NGSolve`](#prerequisite-ngsolve)
  + [Independent `NGSolve` installation](#independent-ngsolve-installation)
    - [Required version of `NGSolve`](#required-version-of-ngsolve)
    - [Installation of `NGSolve`](#installation-of-ngsolve)
  + [Include `NGSolve` installation in `ngsxfem` build.](#include-ngsolve-installation-in-ngsxfem-build)
* [Installation guide for build from source](#installation-guide-for-build-from-source)
  + [Installation steps on Linux](#installation-steps-on-linux)
    - [1a.  Building `ngsxfem` with pre-installed `NGSolve`](#1a-building-ngsxfem-with-pre-installed-ngsolve)
    - [1b.  Building `ngsxfem` and `NGSolve` together](#1b-building-ngsxfem-and-ngsolve-together)
    - [2. Troubleshooting](#2-troubleshooting)
    - [3. Updating `ngsxfem`](#3-updating-ngsxfem)
  + [Installation steps on MacOs](#installation-steps-on-macos)
    - [0.  Building prerequisites on MacOS](#0-building-prerequisites-on-macos)
    - [1.  Building `ngsxfem` with pre-installed `NGSolve`](#1-building-ngsxfem-with-pre-installed-ngsolve)
    - [2.  Updating `ngsxfem`](#2-updating-ngsxfem)
* [Testing the installation](#testing-the-installation)
---


## Prerequisite `NGSolve`
To build `ngsxfem` from source, the corresponding version of `Netgen/NGSolve` is required to be installed. This can either be done in advance (default and recommended option), or as an external dependency. 

### Independent `NGSolve` installation
#### Required version of `NGSolve`
`NGSolve` is pulled as a git submodule. The version of the submodule `ngsolve` in `external_dependencies` is compatible with this version of `ngsxfem`. A necessary version of `NGSolve` can also be found in the first line of [`python/ngs_check.py`](python/ngs_check.py). We always try to have the master branch compatible with the most recent version of the `NGSolve` master branch. If in doubt make sure that you install exactly this version of `NGSolve` before building `ngsxfem`.

#### Installation of `NGSolve`
`NGSolve` can either be [installed from source](https://ngsolve.org/docu/latest/install/install_sources.html) or by installing the latest [pre-built packages (ppa / dmg / ...)](https://ngsolve.org/downloads). Note that the pre-built packages are available in a various older versions. 
Make sure that all environment variables have been set correctly, see [here (for linux)](https://ngsolve.org/docu/latest/install/installlinux.html#finishing-the-installation) and [here (for MacOs)](https://ngsolve.org/docu/latest/install/gettingstarted.html#mac-os-x).

### Include `NGSolve` installation in `ngsxfem` build.
If you have no pre-installed version of `NGSolve` you can include the installation of `NGSolve` into the `ngsxfem` build. In this case `NGSolve` is installed in the version prescribed by the `ngsolve` git submodule. Note that building `NGSolve` will typically take a considerably longer time than `ngsxfem`.

## Installation guide for build from source
### Installation steps on Linux

#### 1a.  Building `ngsxfem` with pre-installed `NGSolve`
Choose a directory where you wish to download the source files and build the library. We shall refer to this location as `BASEDIR`. Here the git repository needs to be cloned.

``` {.shell}
export BASEDIR=`pwd`
git clone https://github.com/ngsxfem/ngsxfem.git src-xfem
```

You then need to create a build directory, configure the build, build and install the build. Here `INSTLOCATION` should be the install directory of `NGSolve`. This depends on the way in which `NGSolve` was installed.

``` {.shell}
mkdir build-xfem
cd build-xfem

cmake \
    -DCMAKE_INSTALL_PREFIX=INSTLOCATION \
    -DBUILD_NGSOLVE=OFF \
    ${BASEDIR}/src-xfem

make
make install
```

You may want to add `-jx` with \'x\' the number of threads you wish to compile with.

#### 1b.  Building `ngsxfem` and `NGSolve` together

If you do not have `Netgen/NGSolve` installed in advance, you can build this as a sub-module. Again, choose a directory where you wish to build and install  everything.

``` {.shell}
export BASEDIR=`pwd`

git clone https://github.com/ngsxfem/ngsxfem.git src-xfem
cd src-xfem
git submodule update --init

cd ${BASEDIR}
mkdir -p ${BASEDIR}/build-xfem ${BASEDIR}/inst

cd ${BASEDIR}/build-xfem
cmake \
    -DCMAKE_INSTALL_PREFIX=${BASEDIR}/inst \
    -DBUILD_NGSOLVE=ON \
    ${BASEDIR}/src-xfem

make
make install
```

Now to start `Netgen` from the command line `${BASEDIR}/inst/bin` has to added to the `PATH`. To run python scripts, the `PYTHONPATH` must be set appropriately

``` {.shell}
export PYTHONPATH=${BASEDIR}/inst/`python3 -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))"`
```

#### 2. Troubleshooting

If you have problems compiling problems or at run time some `NGSolve` symbols are not found, it may be (happened in some configurations) that the `NGSolve` compiler and linker wrapper `ngscxx` and `ngsld` were not used. In this case you may want to add

``` {.shell}
-DCMAKE_CXX_COMPILER=ngscxx -DCMAKE_LINKER=ngsld
```

to the cmake configuration.

#### 3. Updating `ngsxfem`

To update `ngsxfem`, update the source files and build everything
again:

``` {.shell}
cd ${BASEDIR}/src-xfem
git pull

cd ${BASEDIR}/build

make
make install
```

If `NGSolve` was built as a submodule, then after pulling the latest `ngsxfem` sources, also update `NGSolve` by calling `git submodule update --init` in the `src-xfem` directory.

### Installation steps on MacOS

#### 0.  Building prerequisites on MacOs
To build on MacOS you require the Xcode Command Line Tools. These can be installed by calling `xcode-select --install` from within a terminal. Furthermore, CMake must be downloaded and installed. This can be done via [CMake website](https://cmake.org). To use cmake from a terminal, make sure to install the command line tools: Open CMake, in the \"Tools\" menu click on \"How to Install For Command Line Use\" and follow one of the suggested options.

#### 1.  Building `ngsxfem` with pre-installed `NGSolve`

The only difference compared to Linux is that CMake needs to be given the location of the `NGSolve` cmake configuration. This is done by giving the additional flag `-DNGSolve_DIR=INSTLOCATION/Contents/Resources/CMake`. If you have installed `NGSolve` using the dmg file, then `INSTLOCATION` is `/Applications/Netgen.app`. Once `NGSolve` is successfully installed, then `ngsxfem` can be build using the following steps:

``` {.shell}
export BASEDIR=`pwd`
git clone https://github.com/ngsxfem/ngsxfem.git src-xfem

mkdir -p ${BASEDIR}/build-xfem
cd ${BASEDIR}/build-xfem

cmake \
 -DCMAKE_INSTALL_PREFIX=INSTLOCATION \
 -DNGSolve_DIR=INSTLOCATION/Contents/Resources/CMake \
 -DBUILD_NGSOLVE=OFF \
 ${BASEDIR}/src-xfem

make
make install
```

#### 2. Updating `ngsxfem`
See **3. Updating `ngsxfem`** for Linux above.


## Testing the installation

We run tests by default. I you wish to test your self-built binaries, go to the `build-xfem` directory and run `make test` or `ctest`. If you need to see specific tests failing use ctest -V. To run individual tests use ctest -R \<regex\>. E.g. ctest -R cutint to only run cut integration tests. Note that we use `pytest` and `psutil` (with python version \> 3). These can easily be installed using `pip`.


# Docker container
A convenient and reproducible way to set up `ngsxfem` is the usage of [a docker image](https://hub.docker.com/r/ngsxfem/ngsxfem) that we provide here:
<https://hub.docker.com/r/ngsxfem/ngsxfem>.
Installation of `docker` on the common platforms is described [here](https://docs.docker.com/get-docker/). After installation the `Docker daemon` has to be started. This can either be done on boot or manually. In most Linux distributions the command for the latter is either `systemctl start docker` or `service docker start`.

Assuming an installed `docker` with running docker daemon, you can spawn into the `ngsxfem` image with
``` {.shell}
docker run -i -t ngsxfem/ngsxfem:latest /bin/bash
```
## Jupyter through docker
To directly spawn a jupyter server in the docker that you can access from a browser start
``` {.shell}
docker run -p 8888:8888 ngsxfem/ngsxfem-jupyter
```
and open a browser and past in the URL that you obtain in the terminal. You will have the jupyter tutorial files from the docker container available to work with. Note that changes will not be persistent in the image. To work on local files (with persistent changes) mount a local directory to the docker container, e.g.
``` {.shell}
docker run -p 8888:8888 -v ${PWD}:/home/jovyan ngsxfem/ngsxfem-jupyter
```
