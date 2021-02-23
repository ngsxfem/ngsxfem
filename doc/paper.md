---
title: '`ngsxfem`: An add-on to NGSolve for unfitted finite element discretizations'
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
    affiliation: 1
  - name: Janosch Preuß
    orcid: 0000-0002-7384-535X
    affiliation: 1
  - name: Henry von Wahl
    orcid: 0000-0002-0793-1647
    affiliation: 2
affiliations:
 - name: Institut für Numerische und Angewandte Mathematik, Georg-August Universität Göttingen
   index: 1
 - name: Institut für Analysis und Numerik, Otto-von-Guericke Universität, Magdeburg
   index: 2
date: xx xx 2021
bibliography: doc/lit-ngsxfem.bib

---

# Summary
`ngsxfem` is an Add-on library to Netgen/NGSolve which enables the use of unfitted finite element technologies known as XFEM, CutFEM, TraceFEM, Finite Cell,.. . `ngsxfem` is an academic software. Its primary intention is to facilitate the development and validation of new numerical methods.

# Statement of need

Features
* level set described geometries
* element markers
* extended FESpaces
* Cut domain integration
* cut stabilization techniques (ghost penalty)

Distinct features:
* Unfitted space-time (2D and 3D + t) finite elements and quadrature
* isoparametric unfitted FEM
* Multigrid CutMG
* multiple level sets

Applications:
* PDEs on moving domains @Pre18 @Hei20 @Lud20 @LO19 (e.g. two-phase flows, FSI @vWRL20 @vWRFH21) 
* flows on surfaces @JR19a, @JR19, 
* interface problems (Stokes @LPWL16 / scalar @LR16, @Leh16)
* fictitous domains @Leh16a
* model order reduction with geometry 
* particle-laden flows
* shape optimization @Rau18


## Numerical integration on implicitly described (via a level set function) geometries which are not fitted to the mesh
Given a level set function $\phi$ which describes the geometry (e.g. $\Omega = \{ \phi < 0 \}$) a piecewise linear approximation is made.
On simplices (triangles and tetrahedra) this gives a planar intersection on every element which allows for an explicit decomposition into simple geometries.
On these simple (uncut) geometries standard quadrature rules of arbitrary order can be applied which results in quadrature rules for the (approximated) sub-domains where the level set is positive/negative/zero.

![Left: Subdivision strategy for tetrahedra. Right: Integration points on a cut element](doc/graphics/cuttet-quadrature.png){ height=2.5cm align=center }

## Tools to work on an "active mesh" only
In unfitted finite element methods some functions and integrals are only defined on a subset of the mesh. Accordingly finite element spaces and integrals have to be defined only on this active part of the mesh. 
`ngsxfem` offers the tools to mark the corresponding elements and facets and use the marking during assembly and definition of finite element spaces. 
On cut elements one often also uses locally modified finite elements, e.g. by restriction of finite elements on the background mesh.

![unfmesh](doc/graphics/unfittedmesh.jpg){ height=130px} ![xfem](doc/graphics/xfem.jpg){ height=130px}

Left: Active elements with respect to level set function. Right: XFEM basis function.

## Higher order representation of implicit level-set geometries 
To obtain higher order accuracy of integrals, we offer a mesh transformation technique in the spirit of isoparametric finite element methods. 
Thereby the piecewise linear approximation of the level-set (which is only of second order) is mapped onto a higher order accurate approximation of the true geometry.

![Left: Piecewise linear approximation. Right: Higher-order mapped domain](doc/graphics/lsetcurv.jpg){ height=2.5cm align=center} 

## Space-Time Finite Elements for the treatment of moving domain problems
To obtain robust method for partial differential equations on unfitted moving domain we can formulate space-time discretizations. `ngsxfem` provides the necessary tools (so far only in two space dimensions) to define space-time finite element spaces and to integrate on space-time domains. Furthermore, it extends the tools for higher order accurate geometry handling into the space-time setting.

![Left: Sketch of a space-time moving domain (1D + time) Right: sketch of an isoparametrically mapped space-time prism cut by the zero level of a space-time level set function (red).](doc/graphics/spacetime.png){ height=2.5cm align=center}


## Geometries described by multiple level sets
Many geometries do not have -- due to sharp corners or edges -- an implicit description by a *smooth* level set function. Instead multiple level set functions can be used to describe theses geometries. 
To work with these more complicated domains, `ngsxfem` provides tools to work with these geometries as with simple geometries, for example computing the level-set description of the boundary and exterior.  
To enable integration on such domains `ngsxfem` generates quadrature rules with respect to every level-set which cuts a given element. Furthermore, it provides the analogous tools to the single level-set setting, to mark those elements of the mesh which are relevant to a given geometry.

![Left: Elements marked with respect to multiple level sets. Right: Quadrature for multiple cuts.](doc/graphics/mlset.png){ height=3cm align=center}



# `ngsxfem` in the scientific literature
`ngsxfem` has been used in a variety of applications. These include surface problems @JR19a, fluid-structure interaction problems with contact @vWRFH21, reduced order methods @KNBR20 and optimal control problems @AK20. 

As far as we are aware `ngsxfem` has been used in the following scientific literature:

* *Journal Publications and pre-prints*: @vWRFH21, @OQS21, @BR20, @KNBR20, @AK20, @vWRL20, @KBR20, @LR19, @JR19a, @JR19, @LO19, @HL19, @GLR18, @LR17, @Leh17, @LPWL16, @LR16, @Leh16a, @Leh16
* *Accompanying source code* @vWRFH20a, @vWRL20a
* *Thesis* @Hei20, @Lud20, @Hei18, @Rau18, @Pre18,

# Acknowledgements
T. Ludescher, the NGSolve crew!

# References
