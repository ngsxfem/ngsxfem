# Quadrature on cut elements: standard rules and Algoim

`ngsxfem` has two families of quadrature rules for elements (facets, boundary elements, space-time prisms) that are cut by the zero level of one or several level set functions:

* The **standard rules** work with a level set function that is piecewise linear (multilinear on quadrilaterals and hexahedra), given as a `GridFunction` in `H1(mesh, order=1)`. A cut simplex is decomposed into simplices, on which standard rules are applied. Higher order accuracy w.r.t. a smooth geometry is obtained with the isoparametric mapping (`LevelSetMeshAdaptation`).
* The **Algoim-based rules** work with the level set function itself, represented as a polynomial of higher degree on each element. They use the quadrature algorithm of the bundled [Algoim](https://algoim.github.io/) library by R. I. Saye (*J. Comput. Phys.* 448, 110720, 2022). No piecewise linear approximation and no mesh deformation is involved.

Both are used through the same interfaces (`dCut`, `Integrate`, `CutInfo`, ...). Which one is applied is decided for every integral by the option `cut_quadrature` and the type of the level set function.

## The option `cut_quadrature`

```python
from xfem import *

ngsxfemglobals.cut_quadrature = "auto"             # global setting (default)
dx = dCut(lset, NEG, cut_quadrature="algoim")      # for one integral
ci = CutInfo(mesh, lset, cut_quadrature="algoim")  # for the element markers
# dictionary interface
lsetdom = {"levelset": lset, "domain_type": NEG, "cut_quadrature": "algoim"}
```

The values are (as strings or as the enum `CUT_QUADRATURE.AUTO`, `.STANDARD`, `.ALGOIM`):

* `"auto"` (default): the standard rules for a piecewise linear level set `GridFunction`, the Algoim-based rules for all other level set functions and for cases without standard rules.
* `"standard"`: always the standard rules. An exception is raised if they cannot treat the level set function.
* `"algoim"`: the Algoim-based rules wherever they are available, also for piecewise linear level set functions.

A value that is given for an integral or a `CutInfo` takes precedence over the global setting. Element markers and integrals should use the same method.

## Which rule is used?

In the tables, GF stands for `GridFunction` and CF for `CoefficientFunction`. A level set function counts as *P1* if it is a GF in `H1(mesh, order=1)`, in `Discontinuous(H1(mesh, order=1))` or, for space-time, in `ScalarTimeFE(k_t) * H1(mesh, order=1)` (any order in time).

**One level set function.** This covers element integrals, boundary integrals (`vb=BND`), facet integrals (`skeleton=True`), element boundary integrals (`element_boundary=True`) and `CutInfo`.

| level set | `"auto"` | `"standard"` | `"algoim"` |
|:--|:--|:--|:--|
| P1 GF | standard | standard | Algoim |
| P1 GF, prisms | Algoim | exception | Algoim |
| higher order GF | Algoim | exception [1] | Algoim |
| other CF | Algoim | exception [1] | Algoim |

A P1 GF keeps the standard rules also with a mesh deformation, i.e. in the isoparametric approach. The row for prisms refers to element integrals; the triangular and quadrilateral facets of prisms have standard rules.

[1] With `subdivlvl > 0` (element integrals on simplices, no time), `"standard"` uses the subdivision approach: the level set function is interpolated piecewise linearly on a subdivided element.

**Space-time.** This covers `time_order >= 0` and fixed reference times `tref`. The level set function is a GF on a space-time `FESpace` or a CF of the reference time `tref`.

| level set | `"auto"` | `"standard"` | `"algoim"` |
|:--|:--|:--|:--|
| P1 in space | standard | standard | Algoim |
| higher order | Algoim | exception | Algoim |
| CF of `tref` | Algoim | exception | Algoim |

**Multiple level set functions.** This covers a tuple of GFs with a tuple of domain types or a `DomainTypeArray`, and `MultiLevelsetCutInfo`.

| level sets | `"auto"` | `"standard"` | `"algoim"` |
|:--|:--|:--|:--|
| all P1, trig/tet | standard | standard | Algoim |
| all P1, other | Algoim | exception | Algoim |
| higher order | Algoim | exception | Algoim |
| time or `tref` | Algoim | exception | Algoim |
| co-dim. 2 parts | standard [2] | standard [2] | standard [2] |

"Other" are quadrilaterals, hexahedra and prisms. "Higher order" means that at least one level set function is not P1.

[2] Parts of co-dimension two or higher, e.g. `(IF, IF)`, are only supported by the standard rules, i.e. for P1 level set functions on triangles and tetrahedra. They raise an exception otherwise.

The first Algoim-based rule in a run prints one message in the terminal, with the reference to cite. If Algoim has been chosen automatically (and `cut_quadrature` has not been set explicitly), the message says so.

## How the level set function enters the Algoim-based rules

* A level set `GridFunction` is represented exactly, with its polynomial degree on the element.
* Any other `CoefficientFunction` is interpolated on each element by a polynomial of degree `ngsxfemglobals.algoim_lset_order` (default: 4).
* With a mesh deformation, the level set function is evaluated in the mapped points of the deformed element.

Supported are segments, triangles, quadrilaterals, tetrahedra, hexahedra and prisms. The rules are returned in the same form as the standard rules (reference points and weights), so that everything that is built on cut integration rules works unchanged.

## Parameters

* `order` (of `dCut` or `Integrate`): polynomial degree of the integrand up to which the rule is exact on parts that are bounded by planes. The curved geometry is resolved with at least `k+1` Gauss points per direction for a level set polynomial of degree `k`, also for `order=0`.
* `subdivlvl`: number of uniform subdivisions of a cut element before Algoim is applied. This helps on meshes that are coarse compared to the curvature of the interface. Note that this is a different meaning than for the standard rules, see [1] above.
* `ngsxfemglobals.algoim_lset_order`: interpolation degree for level set functions that are not a `GridFunction`.

## Examples

Higher order level set function instead of the isoparametric mapping:

```python
lset_ho = GridFunction(H1(mesh, order=3))
lset_ho.Set(levelset)
ci = CutInfo(mesh, lset_ho)
dx = dCut(lset_ho, NEG, definedonelements=ci.GetElementsOfType(HASNEG))
ds = dCut(lset_ho, IF, definedonelements=ci.GetElementsOfType(IF))
area = Integrate(CoefficientFunction(1) * dx, mesh)
```

Multiple level set functions of higher order:

```python
lsets = tuple(GridFunction(H1(mesh, order=2)) for i in range(2))
lsets[0].Set(x**2 + y**2 - 1)
lsets[1].Set(y - 0.5)
one = CoefficientFunction(1)
half = Integrate(one * dCut(lsets, (NEG, NEG), order=2), mesh)
arc = Integrate(one * dCut(lsets, (IF, NEG), order=2), mesh)
```

Space-time level set function of higher order in space (no space-time mesh deformation):

```python
lset = GridFunction(ScalarTimeFE(2) * H1(mesh, order=2))
lset.Set(levelset)          # levelset depends on t = told + delta_t * tref
dQ = delta_t * dCut(lset, NEG, time_order=4)
dOmnew = dCut(lset, NEG, tref=1)
```

See also `demos/fictdom_algoim.py` and the jupyter tutorials `intlset.ipynb`, `mlset_basic.ipynb` and `spacetime_unfitted.ipynb`.

## Standard rules or Algoim?

Standard rules with the isoparametric mapping:

* The level set function is a P1 `GridFunction`, and the higher order accuracy comes from a mesh deformation. The deformation enters all integrals, and moving domains need a transfer between differently deformed meshes.
* There are few points per cut element and the construction of the rules is cheap.
* Multiple level sets are supported for P1 level set functions on triangles and tetrahedra, including parts of co-dimension two and higher.

Algoim-based rules:

* The level set function is a `GridFunction` of any order or a general `CoefficientFunction`. The mesh is unchanged.
* Prisms are supported in addition to segments, triangles, quadrilaterals, tetrahedra and hexahedra.
* Multiple level sets are supported for any order, on all element types and in space-time, with at most one interface condition per part.
* The rules have many points per cut element and their construction is expensive, see below.

## Cost and limitations of the Algoim-based rules

* The rules have considerably more points and their construction is more expensive. Rough numbers for a sphere with a level set function of order 2: about 800 to 1000 points per cut tetrahedron and about 200 to 450 per cut hexahedron for the volume part, compared to about 30 to 50 points of the standard rules. In 2D a cut triangle has some tens of points. Space-time rules are again considerably more expensive.
* In 3D, level set polynomials of degree 4 and higher can lead to polynomials of a degree beyond 1000 in Algoim's dimension reduction. This happens on elements where the interface has vertical tangents with respect to several coordinate directions, typically on meshes that are coarse compared to the curvature of the interface. An exception reports this case. Remedies are a lower degree of the level set function (for a `CoefficientFunction`: `ngsxfemglobals.algoim_lset_order = 3`), a finer mesh or `subdivlvl > 0`.
* For multiple level sets the rules use a generic perturbation of relative size 1e-9, so that results are accurate up to about 1e-9.
* Not available: parts of co-dimension two and higher of multiple level set domains with level set functions of higher order, and element boundary integrals with time.

## Changes for existing code

* A level set function that is a P1 `GridFunction` is treated exactly as before.
* A higher order `GridFunction` or a `CoefficientFunction` without `subdivlvl` used to be interpolated piecewise linearly on each cut element without notice. It is now integrated with Algoim, which is more accurate, gives different numbers and is more expensive. For the previous behaviour, interpolate the level set function with `InterpolateToP1`.
* A higher order `GridFunction` or a `CoefficientFunction` with `subdivlvl > 0` used the subdivision approach. With `"auto"` it now uses Algoim on a subdivided element. For the previous behaviour, pass `cut_quadrature="standard"`.

## Citation

If you use the Algoim-based quadrature, please also cite

```
R. I. Saye
High-order quadrature on multi-component domains implicitly defined by multivariate polynomials
Journal of Computational Physics, 448, 110720, 2022,
https://doi.org/10.1016/j.jcp.2021.110720
```
