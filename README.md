# About `ngsxfem`

`ngsxfem` is an add-on library to the finite element package [Netgen/NGSolve](https://ngsolve.org) which enables the use of unfitted finite element technologies known as XFEM, CutFEM, TraceFEM, Finite Cell, ... . `ngsxfem` is an academic software. Its primary intention is to facilitate the development and validation of new numerical methods.

# The features of `ngsxfem`

The main features of `ngsxfem` are:

* Tools to work on an a subset of the trianguation, the \"active mesh\" only

* Numerical integration on geometries that are (implicitly) described through level set functions.

* Higher order representation of level set geometries

* Space-Time Finite Elements for the treatment of moving domain problems

* All these features combined with the usual flexibility and power of [NGSolve](https://ngsolve.org).

`ngsxfem` has been used in a variety of applications. In [`doc/paper.md`](doc/paper.md) more details on the features and references to applications where `ngsxfem` is used are given.

Not all features of `ngsxfem` and `NGSolve` can directly be combined. Here is an overview of `ngsxfem` and `NGSolve` features and if they can directly be combined:
| Features ⇲| `CFE` | `XFE` | `DGF` | `Iso` | `MLS` | `STF` | `GhP` | `Hex` | `Tet` | `MPI` |
|-|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| `CFE`: CutFEM formulation | / | / | yes | yes | yes | no | yes | yes | yes | yes |
| `XFE`: XFEM formulation | / | / | yes | yes | no | no | yes | yes | yes | yes |
| `DGF`: Discontinuous Galerkin form. | yes | yes | / | yes | no | no | yes | no | yes | no |
| `Iso`: isoparametric mesh deformation | yes | yes | yes | / | no | yes | yes | yes | yes | yes |
| `MLS`: multiple level set | yes | no | no | no | / | no | yes | no | yes | yes |
| `STF`: space-time finite elements | yes | no | no | yes | no | / | yes | yes | yes | yes |
| `GhP`: Ghost penalty | yes | yes | yes | yes | yes | yes | / | yes | yes | no |
| `Hex`: hyperrectangles (quads/hexes) | yes | yes | no? | yes | no | yes | yes | / | / | yes |
| `Tet`: simplices (trigs./tets) | yes | yes | yes | yes | yes | yes | yes | / | / | yes |
| `MPI`: MPI | yes | yes | no | yes | yes | yes | no | yes | yes | / |

Some of the *no*s are work in progress (e.g. `MLS`&`STF`) and some have not been considered so far (e.g. `DGF`&`STF`). If you need a certain combination to work, please contact us and we will see what we can do. 

# Examples and Documentation

We provide two main sources to learn how to use `ngsxfem`:
 * At <https://github.com/ngsxfem/ngsxfem-jupyter> you can find tutorial-style jupyter notebooks for ngsxfem. These explain the core functionalities and usage of the tools provided by `ngsxfem`. You can run those tutorials interactively (without the need of a local installation) through binder. 
 * in the [`demos`](./demos)-directory we provide several examples that demonstrate the usage of `ngsxfem` features. See [`demos/README.md`](demos/README.md) for details.

 # Installation
 We provide installation instructions for installation through `pip` and installation from source in [`INSTALLATION.md`](INSTALLATION.md). Further, [a docker image](https://hub.docker.com/r/schruste/ngsxfem) is available that can be used to run `ngsxfem` through docker.

 # List of contributing authors (with major contributions)

-   Christoph Lehrenfeld (main author)
-   Fabian Heimann (cutIntegration, space-time)
-   Thomas Ludescher (multigrid)
-   Janosch Preuss (space-time)
-   Henry von Wahl (multiple levelsets)
