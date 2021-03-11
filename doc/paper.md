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

# Scope and purpose
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

`ngsxfem` has already been used to develop and investigate the performance of different unfitted discretizations in several application fields:

* Scalar fictitious domain problems have been investigated in @Leh17.
* Scalar and Stokes interface problems as in two-phase flow problems have been investigated in @LR16, @Leh16, @Leh16a, @LR19, @LPWL16, @OQS21 and @Lud20. In the latter publication Multigrid preconditioners have been developed. 
* PDEs on moving domains have been considered and discretizations have been implemented in `ngsxfem` in several publications. In @Pre18, @Hei20 space-time discretizations for scalar parabolic model problems are considered, whereas Eulerian time stepping methods for scalar and unsteady Stokes problems have been considered in @LO19 and @vWRL20, see also the reproduction data in @vWRL20a. 
* Fluid-structure interaction problems have been investigated in @vWR21 and @vWRFH21, see also the reproduction data in @vWRFH20a.
* Surface PDEs have also been considered ranging from scalar ones in @GLR18 and @Hei18 over Vector-Laplacians in @JR19 and @R20 to flows on smooth surfaces in @JR19a and @BR20.
* Shape optimization for geometries that are described by level set functions have been considered in @Rau18.
* Model order reduction and optimal control with geometry have been considered in @AK20 and @KBR20 .

## Feature details
In the previous section we briefly explained or mentioned some tools of `ngsxfem`. In this section we elaborate on a few non-trivial aspects of them.

### Tools to work on an "active mesh" only
In unfitted finite element methods some functions and integrals are only defined on a subset of the mesh. Accordingly finite element spaces and integrals have to be defined only on this active part of the mesh. 
`ngsxfem` offers the tools to mark the corresponding elements and facets and use the marking during assembly and definition of finite element spaces. 
On cut elements one often also uses locally modified finite elements, e.g. by restriction of finite elements on the background mesh.

![Left: Active mesh. Right: Basis function restricted to a subdomain.](graphics/xfes.png){ height=3cm align=center}

### Numerical integration on unfitted geometries described by one level set function
Given a level set function $\phi$ which describes the geometry (e.g. $\Omega = \{ \phi < 0 \}$), a piecewise linear (or bilinear on hyperrectangles, cf. @HL19) approximation is made. How to obtain a higher order reconstruction from this basis approximation is discussed [below](#higher-order-representation-of-implicit-level-set-geometries).

On simplices (triangles and tetrahedra) this gives a planar intersection on every element which allows for an explicit decomposition into simple geometries.
On these simple (uncut) geometries standard quadrature rules of arbitrary order can be applied which results in quadrature rules for the (approximated) sub-domains where the level set is positive/negative/zero.

![Left: Subdivision strategy for cut tetrahedron. Right: Integration points on a cut element.](graphics/cuttet-quadrature.png){ height=2.5cm align=center }



### Geometries described by multiple level sets
While one level set function may be sufficient for the approximation of a smooth geometry, many geometries do not have -- due to sharp corners or edges -- an implicit description by only one *smooth* level set function. In these cases often multiple level set functions can be used to describe theses geometries. 
To work with these more complicated domains, `ngsxfem` provides tools to work with these geometries as with simple geometries, for example computing the level set description of the boundary and exterior.  To enable integration on such domains `ngsxfem` generates quadrature rules with respect to every level set which cuts a given element. Furthermore, it provides the analogous tools to the single level-set setting, to mark those elements of the mesh which are relevant to a given geometry.

![Left: Elements marked w.r.t. multiple level sets. Right: Quadrature for multiple cuts.](graphics/mlset.png){ height=3cm align=center}


### Higher order representation of implicit level-set geometries 
To obtain higher order accuracy in the handling of the geometry, `ngsxfem` offers a mesh transformation technique in the spirit of isoparametric finite element methods. 
Thereby the piecewise linear approximation of the level set (which is only of second order) is mapped onto a higher order accurate approximation of the true geometry.

![Left: Piecewise linear approximation. Right: Higher-order mapped domain](graphics/lsetcurv.jpg){ height=2.5cm align=center} 

### Space-Time Finite Elements for the treatment of moving domain problems
To obtain robust methods for partial differential equations on unfitted moving domain we can formulate space-time discretizations. `ngsxfem` provides the necessary tools to define space-time finite element spaces and to integrate on space-time domains. Furthermore, it extends the tools for higher order accurate geometry handling into the space-time setting.

![Left: Sketch of a space-time moving domain (1D + time) Right: sketch of an isoparametrically mapped space-time prism cut by the zero level of a linear-in-space space-time level set function on the reference configuration (red).](graphics/spacetime.png){ height=2.5cm align=center}

# Acknowledgements
The authors acknowledge the support by the `NGSolve` crew, especially Matthias Hochsteger and Christopher Lackner for keeping the build system compatible with `NGSolve`. The authors also want to thank Thomas Ludescher for his developments on MultiGrid for unfitted FEM that he contributed to the project. Further, part of the implementation of the numerical integration routines has been developed within the project “LE 3726/1-1” funded by the Deutsche Forschungsgemeinschaft (DFG, German Science Foundation). Henry von Wahl has been funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - 314838170, GRK 2297 MathCoRe.

# References
