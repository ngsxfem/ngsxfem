[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/ngsxfem/ngsxfem/pypi.yml?branch=master&label=build%20and%20pypi&logo=github)](https://github.com/ngsxfem/ngsxfem/actions/workflows/pypi.yml) 
 [![GitHub release (latest by date)](https://img.shields.io/github/v/release/ngsxfem/ngsxfem?label=latest%20release&logo=github)](https://github.com/ngsxfem/ngsxfem) 
[![jossstatus](https://joss.theoj.org/papers/9fda1eadfc58af64b89dc7f27043f4cb/status.svg)](https://joss.theoj.org/papers/9fda1eadfc58af64b89dc7f27043f4cb)

 

[![PyPI](https://img.shields.io/pypi/v/xfem?color=blue&label=latest%20PyPI%20version&logo=pypi&logoColor=white)](https://pypi.org/project/xfem/)
[![PyPI - Implementation](https://img.shields.io/pypi/implementation/xfem?logo=pypi&logoColor=white)](https://pypi.org/project/xfem/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/xfem?color=green&label=PyPI%20downloads&logo=pypi&logoColor=white)](https://pypi.org/project/xfem/)

[![Docker Pulls](https://img.shields.io/docker/pulls/ngsxfem/ngsxfem?logo=docker&logoColor=white)](https://hub.docker.com/repository/docker/ngsxfem/ngsxfem) 
[![badge](https://img.shields.io/badge/launch-jupyter--tutorials-579ACA.svg?logo=data:image/png;base64,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)](https://mybinder.org/v2/gh/ngsxfem/ngsxfem-jupyter/HEAD?filepath=tutorials.ipynb)

[![lite-badge](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://ngsuite.pages.gwdg.de/ngsxfem-jupyter/lab/index.html)

### Documentation

A *preliminary* collection of documentation (demos, tutorials, mini API) can be found [here](https://ngsuite.pages.gwdg.de/ngsxfem/).
A simple doxygen collection of `NGSolve`, `ngsxfem` and `ngstrefftz` C++ API can be found [here](https://ngsuite.pages.gwdg.de/ngs-doxygen/).

# About `ngsxfem`

`ngsxfem` is an add-on library to the finite element package [Netgen/NGSolve](https://ngsolve.org) which enables the use of unfitted finite element technologies known as XFEM, CutFEM, TraceFEM, Finite Cell, ... . `ngsxfem` is an academic software. Its primary intention is to facilitate the development and validation of new numerical methods for partial differential equations.

# The features of `ngsxfem`

The main features of `ngsxfem` are:

* Tools to work on a subset of the triangulation, the \"active mesh\" only
* Numerical integration on geometries that are (implicitly) described through level set functions.
* Higher order representation of level set geometries
* Space-Time Finite Elements for the treatment of moving domain problems
* All these features combined with the usual flexibility and power of [NGSolve](https://ngsolve.org).

`ngsxfem` has been used in a variety of applications. In the document [doc/feature-details.md](doc/feature-details.md) (see also [compiled pdf](https://nightly.link/ngsxfem/ngsxfem/workflows/extras-workflow/master/doc-features.zip) ) more details on the features is given and in [doc/literature.md](doc/literature.md) (see also [literature](https://nightly.link/ngsxfem/ngsxfem/workflows/extras-workflow/master/doc-literature.zip) ) an overview of the scientific literature where `ngsxfem` is used is provided.

Not all features of `ngsxfem` and `NGSolve` can directly be combined. Here is an overview of `ngsxfem` and `NGSolve` features and if they can directly be combined:

| Features ⇲| `CFE` | `XFE` | `DG` | `Iso` | `MLS` | `ST` | `Gh` | `Ag` | `Hex` | `Tet` | `MPI` |
|-|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| `CFE`: CutFEM form. | / | / | y | y | y | y | y | y | y | y | y |
| `XFE`: XFEM formulation | / | / | y | y | n | n | y | n | y | y | y |
| `DG`: Discont. Galerkin  | y | y | / | y | n | y| y | y| y| y | n |
| `Iso`: isoparametric map | y | y | y | / | n | y | y | y| y | y | y |
| `MLS`: multiple level set | y | n | n | n | / | n | y | y| n | y | y |
| `ST`: space-time FEM | y | n | y | y | n | / | y | n | y | y | y |
| `Gh`: Ghost penalty | y | y | y | y | y | y | / | / |  y | y | n |
| `Ag`: Agg. FEM | y | n | y | y | y | n | / | / | y | y | n |
| `Hex`: quads / hexes | y | y | y | y | n | y | y | y | / | / | y |
| `Tet`: trigs./tets | y | y | y | y | y | y | y | y | / | / | y |
| `MPI`: MPI | y | y | n | y | y | y | n | n | y | y | / |

Some of the *no*s are work in progress and some have not been considered so far. If you need a certain combination to work, please contact us and we will see what we can do. 

# Examples and Documentation

We provide two main sources with which to learn how to use `ngsxfem`:
 * At <https://github.com/ngsxfem/ngsxfem-jupyter> you can find tutorial-style jupyter notebooks for `ngsxfem`. These explain the core functionalities and usage of the tools provided by `ngsxfem`. You can run those tutorials interactively (without the need of a local installation) through binder. 
 * In the [demos](./demos)-directory we provide several examples that demonstrate the usage of `ngsxfem` features. See [demos/README.md](demos/README.md) for details.

# Installation
 We provide installation instructions for building/installing through `pip` and building/installing from sources in [INSTALLATION.md](INSTALLATION.md). Further, [a docker image](https://hub.docker.com/r/ngsxfem/ngsxfem) is available which can be used to run `ngsxfem` through docker.

# List of contributing authors

Major contributions:

* Christoph Lehrenfeld (main author)
* Fabian Heimann (cut integration, space-time, AggFEM)
* Henry von Wahl (multiple levelsets, mac support, AggFEM)
* Janosch Preuss (space-time)
* Thomas Ludescher (multigrid)
* Paul Stocker (CI, docu, builds, AggFEM)

Additional contributions:

* Pedro Costa Klein (CI, docu)
* Maximilian Zienecker (SIMD, CI)

# Community guidelines
If you observe any problems with the software / examples / documentation / installation or want to contribute, you can get in touch with us through either:
* the [github issue tracker](https://github.com/ngsxfem/ngsxfem/issues) 
* the [NGSolve user forum](https://forum.ngsolve.org/tag/ngsxfem).

# Citing

If you use `ngsxfem` for academic work, please consider citing our publication:

```
C. Lehrenfeld, F. Heimann, J. Preuß and H. von Wahl
ngsxfem: Add-on to NGSolve for geometrically unfitted finite element discretizations
Journal of Open Source Software, 6(64), 3237,
https://doi.org/10.21105/joss.03237
```

