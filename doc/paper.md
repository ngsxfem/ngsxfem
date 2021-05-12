---
title: '`ngsxfem`: Add-on to NGSolve for geometrically unfitted finite element discretizations'
tags:
  - finite elements
  - unfitted finite elements
  - level set geometry
  - Cut FEM
  - eXtended FEM
  - Cut-Cell methods
authors:
  - name: Christoph Lehrenfeld^[Corresponding Author]
    orcid: 0000-0003-0170-8468
    affiliation: 1
  - name: Fabian Heimann
    orcid: 0000-0002-8969-6504
    affiliation: 1
  - name: Janosch Preuß
    orcid: 0000-0002-7384-535X
    affiliation: 1
  - name: Henry von Wahl
    orcid: 0000-0002-0793-1647
    affiliation: 2
affiliations:
 - name: Institute of Numerical and Applied Mathematics, Georg-August Universität Göttingen
   index: 1
 - name: Institute for Analysis and Numerics, Otto-von-Guericke Universität, Magdeburg
   index: 2
date: xx xx 2021
bibliography: lit-ngsxfem.bib
---

# Summary
`ngsxfem` is an add-on library to [`Netgen/NGSolve`](www.ngsolve.org), a general purpose, high performance finite element library for the numerical solution of partial differential equations. The add-on enables the use of geometrically unfitted finite element technologies known under different labels, e.g. *XFEM*, *CutFEM*, *TraceFEM*, *Finite Cell*, *fictitious domain method* or *Cut-Cell methods*, etc.. Both, `Netgen/NGSolve` and `ngsxfem` are written in C++ with a rich python interface through which it is typically used. `ngsxfem` is an academic software. Its primary intention is to facilitate the development and validation of new numerical methods.

# Statement of need
Typically, in the finite element method for the discretization of PDEs, the geometry under consideration is parametrized by the computational mesh yielding *geometrically fitted* finite element methods. The generation and adaptation of geometrically fitted computational meshes can be a burden on simulation methods, e.g. if the geometries are complex or especially if they are evolving in time. To be more flexible w.r.t. the geometry, *geometrically unfitted* finite element methods can be considered which break the direct link between the geometry parametrization and the computational mesh. Instead, a separate description of the geometry, e.g. through a *level set function* is used. `ngsxfem` aims at providing the necessary tools to robustly work in a *geometrically unfitted* setting where the geometry is described by one (or multiple) *level set function(s)*. The essential tools extending standard finite element codes for the *geometrically unfitted* setting are:

* formulation of geometrically unfitted geometry through level set function(s)
* classification of elements in the computational mesh according to the unfitted geometry
* finite element spaces that consider the cut information 
* facilities to provide robust numerical integration on cut elements 
* stabilization techniques to deal with arbitrary bad cuts (e.g. "ghost penalty")

First of all `ngsxfem` provides these tools for `Netgen/NGSolve`. For other finite element frameworks similar libraries exists, e.g. `libcutfem` or `multimesh` for the `FEniCS` project, cf. @fenics and `dune-udg`, cf. @dune:udg, for `dune`, cf. @dune (more precisely `dune-pdelab`, cf. @dune:pdelab ). Let us also mention the `FEMPAR` finite element package which directly handles unfitted geometries, cf.  @fempar. In addition, `ngsxfem` has three advanced features beyond that:

* higher order handling of curved level set geometries using isoparametric unfitted FEM, cf. @Leh16, @LR16, @Leh17
* space-time finite elements and quadrature for unfitted space-time finite element discretizations of PDEs on moving domains, cf. @Pre18, @Hei20
* the so-called direct version of the ghost penalty stabilization as introduced in @Pre18

A more detailed overview of the key features provided by `ngsxfem` can be found in the repository under `doc/feature-details.md` Furthermore, a detailed overview of the scientific literature which has utilised `ngsxfem` is given in the repository under `doc/literature.md`


# Acknowledgements
The authors acknowledge the support by the `NGSolve` crew, especially Matthias Hochsteger and Christopher Lackner for keeping the build system compatible with `NGSolve`. The authors also want to thank Thomas Ludescher for his developments on MultiGrid for unfitted FEM that he contributed to the project. Further, part of the implementation of the numerical integration routines has been developed within the project “LE 3726/1-1” funded by the Deutsche Forschungsgemeinschaft (DFG, German Science Foundation). Henry von Wahl has been funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - 314838170, GRK 2297 MathCoRe.

# References
