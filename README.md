# About `ngsxfem`

`ngsxfem` is an add-on library to the finite element package [Netgen/NGSolve](https://ngsolve.org) which enables the use of unfitted finite element technologies known as XFEM, CutFEM, TraceFEM, Finite Cell, ... . `ngsxfem` is an academic software. Its primary intention is to facilitate the development and validation of new numerical methods.

**TODO:**Add 'manual' Table of content here (ignore  'About ngsxfem' )
# The main features

The main features of this tool are listed below. In [*doc/paper.md*]{.spurious-link target="doc/paper.md"} more details and references are given.

* Numerical integration on implicitly described geometries. More specifically geometries that are describe by one or multiple level set functions and hence are not fitted to the computational mesh.

* Tools to work on an \"active mesh\" only

* Higher order representation of implicit level-set geometries

* Space-Time Finite Elements for the treatment of moving domain problems

* Geometries described by multiple level sets

### feature matrix

Not all features combine directly. Here is an overview of features and which work or may not work together:

| Features | `CFE` | `XFE` | `DGF` | `Iso` | `MLS` | `STF` | `GhP` | `Hex` | `Tet` | `MPI` |
|-|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| `CFE`: CutFEM formulation | / | / | yes | yes | yes | no | yes | yes | yes | yes |
| `XFE`: XFEM formulation | / | / | yes | yes | no? | no | yes | yes | yes | yes |
| `DGF`: Discontinuous Galerkin form. | yes | yes | / | yes | no | no? | yes | no? | yes | no* |
| `Iso`: isoparametric mesh deformation | yes | yes | yes | / | no | yes | yes | yes | yes | yes |
| `MLS`: multiple level set | yes | no? | no | no | / | no | yes | no | yes | yes |
| `STF`: space-time finite elements | yes | no | no? | yes | no | / | yes | yes | yes | yes |
| `GhP`: Ghost penalty | yes | yes | yes | yes | yes | yes | / | yes | yes | no |
| `Hex`: hyperrectangles (quads/hexes) | yes | yes | no? | yes | no | yes | yes | / | / | yes |
| `Tet`: simplices (trigs./tets) | yes | yes | yes | yes | yes | yes | yes | / | / | yes |
| `MPI`: MPI | yes | yes | no* | yes | yes | yes | no | yes | yes | / |


Some of the nos are work in progress (), some have not been considered so far and some are nos by design. 

`ngsxfem` has been used in a variety of applications see also [paper.md](doc/paper.md).

# Installation

## Docker container

``` {.shell}
docker run -i -t schruste/ngsxfem:latest /bin/bash
```

## Installation through `pip` (`pip3`)

**TODO**
### Releases
* Linux
* Mac
* Windows is not supported a.t.m. (try WSL)

### dev-branch (`master`)
* **TODO** pip command: Linux / Mac (make sure that a from-source build is chosen)

## Building from source

To build `ngsxfem` from source, the corresponding version of `Netgen/NGSolve` is required to be installed. This can either be done in advance (default option), or as an external dependency. `ngsolve` is pulled as a submodule. The version of the submodule is compatible with this version of `ngsxfem`. If in doubt make sure that you install exactly this version of `NGSolve` before building `ngsxfem`.

### Linux

#### 1.  Prerequisites
There are no additional dependencies that come on top of `NGSolve`, see [www.ngsolve.org](https://ngsolve.org/docu/latest/install/installlinux.html).

#### 2.  Building `ngsxfem` with pre-installed `NGSolve`
Make sure that the installed version of `NGSolve` is compatible with the current `ngsxfem` release. If you are building the latest release of `ngsxfem`, then the latest release of `NGSolve` should work.

Choose a directory where you wish to download the source files and build the library. We shall refer to this location as `BASEDIR`. Here the git repository need to be cloned.

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

#### 3.  Building the NGS-Suite and `ngsxfem` together

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

#### 4.  Fix of potential issues

If you have compiling problems or at run time some NGSolve symbols are not found, it may be (happened in some configurations) that the NGSolve compiler and linker wrapper `ngscxx` and `ngsld` were not used. In this case you may add

``` {.shell}
-DCMAKE_CXX_COMPILER=ngscxx -DCMAKE_LINKER=ngsld
```

to the cmake configuration.

#### 5.  Updating `ngsxfem`

To update `ngsxfem`, update the source files and build everything
again:

``` {.shell}
cd ${BASEDIR}/src-xfem
git pull

cd ${BASEDIR}/build

make
make install
```

If `NGSolve` was built as a submodule, then after pulling the latest `ngsxfem` sources, also update NGSolve by calling `git submodule update --init` in the `src-xfem` directory.

### MacOS

To build on MacOS you require the Xcode Command Line Tools. These can be installed by calling `xcode-select --install` from within a terminal. Furthermore, CMake must be downloaded and installed. This can be done via [CMake website](https://cmake.org). To use cmake from a terminal, make sure to install the command line tools: Open CMake, in the \"Tools\" menu click on \"How to Install For Command Line Use\" and follow one of the suggested options.

We recommend that you install `NGSolve` this can either be done [from source](https://ngsolve.org/docu/latest/install/installmacnative.html) or by installing the latest [pre-built dmg](https://ngsolve.org/downloads). Make sure that all environment variables have been [set correctly](https://ngsolve.org/docu/latest/install/gettingstarted.html#mac-os-x).

The only difference compared to linux is that CMake needs to be given the location of the NGSolve cmake configuration. This is done by giving the additional flag `-DNGSolve_DIR=INSTLOCATION/Contents/Resources/CMake`. If you have installed NGSolve using the dmg file, then `INSTLOCATION` is `/Applications/Netgen.app`. Once NGSolve is successfully installed, then `ngsxfem` can be build using the following steps:

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

## Testing the installation

We run test by default. I you wish to test your self-built binaries, go to the `build-xfem` directory and run `make test` or `ctest`. If you need to see specific tests failing use ctest -V. To run individual tests use ctest -R \<regex\>. E.g. ctest -R cutint to only run cut integration tests. Note that we use `pytest` and `psutil` (with python version \> 3). These can easily be installed using `pip`.

# Examples and Documentation

We provide two main sources to learn how to use `ngsxfem`:
 * At <https://github.com/ngsxfem/ngsxfem-jupyter> you can find tutorial-style jupyter notebooks for ngsxfem. These explain the core functionalities and usage of the tools provided by `ngsxfem`. You can run those tutorials interactively (without the need of a local installation) through binder. 
 * in the [`demos`](./demos)-directory we provide several examples that demonstrate the usage of `ngsxfem` features:

| filename | description | 
|-|-|
| [`demos/lsetgeoms.py`](demos/lsetgeoms.py) | Examples of smooth domains described by level set functions and their geometrical approximation. |
| [`demos/unf_interf_prob.py`](demos/unf_interf_prob.py) | Stationary scalar unfitted interface problem with an unfitted isoparametrics FEM discretization and Nitsche. |
| [`demos/fictdom.py`](demos/fictdom.py) | Poisson problem on an unfitted (smooth) domain described by a level set function, discretized with CutFEM, Ghost Penalties and Nitsche. |
| [`demos/fictdom_mlset.py`](demos/fictdom_mlset.py) | Poisson problem on an unfitted domain described by multiple level sets, discretized with CutFEM and Nitsche. |
| [`demos/stokescutfem.py`](demos/stokescutfem.py) | Stationary Stokes interface problem with an unfitted isoparametric Taylor-Hood-Nitsche discretization |
| [`demos/tracefem.py`](demos/tracefem.py) | stationary surface PDE problem with a TraceFEM discretization |
| [`demos/fictdom_dg.py`](demos/fictdom_dg.py) |   Poisson problem on an unfitted (smooth) domain described by a level set function, discretized with a Cut Discontinuous Galerkin method and Nitsche.|
| [`demos/moving_domain.py`](demos/moving_domain.py) | Parabolic PDE on a moving domain solved by an Eulerian time stepping method using CutFEM and Ghost penalty stabilization. |

# List of contributing authors (with major contributions)

-   Christoph Lehrenfeld (main author)
-   Fabian Heimann (cutIntegration, space-time)
-   Thomas Ludescher (multigrid)
-   Janosch Preuss (space-time)
-   Henry von Wahl (multiple levelsets)
