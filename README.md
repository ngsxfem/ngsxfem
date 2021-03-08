[![Build, test (and publish release) of ngsxfem](https://github.com/ngsxfem/ngsxfem/actions/workflows/build-and-test-ngsxfem.yml/badge.svg)](https://github.com/ngsxfem/ngsxfem/actions/workflows/build-and-test-ngsxfem.yml) [![Compile JOSS Paper](https://github.com/ngsxfem/ngsxfem/actions/workflows/paper-workflow.yml/badge.svg)](https://github.com/ngsxfem/ngsxfem/actions/workflows/paper-workflow.yml) [![GitHub commit activity](https://img.shields.io/github/commit-activity/y/ngsxfem/ngsxfem)](https://github.com/ngsxfem/ngsxfem) 
[![GitHub last commit](https://img.shields.io/github/last-commit/ngsxfem/ngsxfem)](https://github.com/ngsxfem/ngsxfem) 
[![LGPL](https://img.shields.io/github/license/ngsxfem/ngsxfem)](https://github.com/ngsxfem/ngsxfem/blob/release/LICENSE)


[![PyPI](https://img.shields.io/pypi/v/xfem?color=blue&label=latest%20PyPI%20version)](https://pypi.org/project/xfem/)
[![PyPI - Wheel](https://img.shields.io/pypi/wheel/xfem)](https://pypi.org/project/xfem/)
[![PyPI - Implementation](https://img.shields.io/pypi/implementation/xfem)](https://pypi.org/project/xfem/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/xfem?color=green&label=PyPI%20downloads)](https://pypi.org/project/xfem/)

[![docker build](https://img.shields.io/docker/cloud/build/ngsxfem/ngsxfem)](https://hub.docker.com/repository/docker/ngsxfem/ngsxfem) 
[![Docker Pulls](https://img.shields.io/docker/pulls/ngsxfem/ngsxfem)](https://hub.docker.com/repository/docker/ngsxfem/ngsxfem) 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ngsxfem/ngsxfem-jupyter/HEAD?filepath=tutorials.ipynb)

# About `ngsxfem`

`ngsxfem` is an add-on library to the finite element package [`Netgen`/`NGSolve`](https://ngsolve.org) which enables the use of unfitted finite element technologies known as XFEM, CutFEM, TraceFEM, Finite Cell, ... . `ngsxfem` is an academic software. Its primary intention is to facilitate the development and validation of new numerical methods for partial differential equations.

# The features of `ngsxfem`

The main features of `ngsxfem` are:

* Tools to work on an a subset of the triangulation, the \"active mesh\" only

* Numerical integration on geometries that are (implicitly) described through level set functions.

* Higher order representation of level set geometries

* Space-Time Finite Elements for the treatment of moving domain problems

* All these features combined with the usual flexibility and power of [NGSolve](https://ngsolve.org).

`ngsxfem` has been used in a variety of applications. In [`doc/paper.md`](doc/paper.md) more details on the features and references to applications where `ngsxfem` is used are given.

Not all features of `ngsxfem` and `NGSolve` can directly be combined. Here is an overview of `ngsxfem` and `NGSolve` features and if they can directly be combined:
| Features ⇲| `CFE` | `XFE` | `DGF` | `Iso` | `MLS` | `STF` | `GhP` | `Hex` | `Tet` | `MPI` |
|-|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| `CFE`: CutFEM form. | / | / | yes | yes | yes | no | yes | yes | yes | yes |
| `XFE`: XFEM formulation | / | / | yes | yes | no | no | yes | yes | yes | yes |
| `DGF`: Discont. Galerkin  | yes | yes | / | yes | no | no | yes | no | yes | no |
| `Iso`: isoparametric map | yes | yes | yes | / | no | yes | yes | yes | yes | yes |
| `MLS`: multiple level set | yes | no | no | no | / | no | yes | no | yes | yes |
| `STF`: space-time FEM | yes | no | no | yes | no | / | yes | yes | yes | yes |
| `GhP`: Ghost penalty | yes | yes | yes | yes | yes | yes | / | yes | yes | no |
| `Hex`: quads / hexes | yes | yes | no | yes | no | yes | yes | / | / | yes |
| `Tet`: trigs./tets | yes | yes | yes | yes | yes | yes | yes | / | / | yes |
| `MPI`: MPI | yes | yes | no | yes | yes | yes | no | yes | yes | / |

Some of the *no*s are work in progress (e.g. `MLS`&`STF`) and some have not been considered so far (e.g. `DGF`&`STF`). If you need a certain combination to work, please contact us and we will see what we can do. 

# Examples and Documentation

We provide two main sources with which to learn how to use `ngsxfem`:
 * At <https://github.com/ngsxfem/ngsxfem-jupyter> you can find tutorial-style jupyter notebooks for `ngsxfem`. These explain the core functionalities and usage of the tools provided by `ngsxfem`. You can run those tutorials interactively (without the need of a local installation) through binder. 
 * in the [`demos`](./demos)-directory we provide several examples that demonstrate the usage of `ngsxfem` features. See [`demos/README.md`](demos/README.md) for details.

 # Installation
 We provide installation instructions for building/installing through `pip` and building/installing from sources in [`INSTALLATION.md`](INSTALLATION.md). Further, [a docker image](https://hub.docker.com/r/schruste/ngsxfem) is available which can be used to run `ngsxfem` through docker.

# List of contributing authors (with major contributions)

-   Christoph Lehrenfeld (main author)
-   Fabian Heimann (cut integration, space-time)
-   Thomas Ludescher (multigrid)
-   Janosch Preuss (space-time)
-   Henry von Wahl (multiple levelsets)

# Community guidelines
If you observe any problems with the software / examples / documentation / installation or want to contribute, you can get in touch with us through either:
 * the [github issue tracker](https://github.com/ngsxfem/ngsxfem/issues) 
 * the [`NGSolve` user forum](https://ngsolve.org/forum/index).
