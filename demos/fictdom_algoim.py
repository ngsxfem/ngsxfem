"""
In this example we solve the unfitted Poisson problem of fictdom.py, but
the higher order geometry approximation is realized with the quadrature of
the Algoim library (R. I. Saye, https://algoim.github.io/) instead of the
isoparametric approach.

Domain:
-------
The domain is [-1,1]^2 while the physical domain is a ring described by a
(polynomial) level set function. In the discretisation, the level set
function is interpolated into a finite element space of order k (no
piecewise linear approximation and no mesh deformation). The cut
integrals are computed w.r.t. this higher order level set function with
Algoim's quadrature (ngsxfemglobals.cut_quadrature, see below).

Compared to the standard approach:
* No mesh deformation (LevelSetMeshAdaptation) is needed, the level set
  function is simply interpolated (or given as CoefficientFunction).
* The same code works on quadrilateral (and hexahedral) meshes.
* For a piecewise linear level set function (lset_order = 1), the standard
  rules (straight cuts) are used, so that the geometry error and the
  discretisation error converge only with second order. With a level set
  function of order k, the L2 error converges with order k+1.

PDE problem:
------------
    - (u_xx + u_yy) = f in Omega (where levelset is negative)
                  u = 0 on dOmega (where levelset is zero)

The r.h.s. term f is chosen according to a manufactured solution.

Discretisation:
---------------
* Background finite element space restricted to active domain (CutFEM)

* Nitsche formulation to impose boundary conditions, see. e.g. [1]

* Ghost penalty stabilization to deal with bad cuts (version as in [2])

Implementational aspects:
-------------------------
* Quadrature on the higher order level set domain with Algoim [3]. The
  quadrature method is selected with ngsxfemglobals.cut_quadrature:
  "standard", "algoim" or "auto" (Algoim for higher order level sets
  without mesh deformation). It can also be chosen for single integrals,
  e.g. dCut(lset, NEG, cut_quadrature="algoim").

* A (sparse) direct solver is applied to solve the arising linear systems.

References:
-----------
All concepts that are used here are explained in the jupyter tutorials
`basics.ipynb`, `intlset.ipynb` and `cutfem.ipynb`.

Literature:
-----------
[1] E. Burman, P. Hansbo, Fictitious domain finite element methods using
    cut elements: II. A stabilized Nitsche method, Appl. Num. Math.
    62(4):328-341, 2012.
[2] J. Preuß, Higher order unfitted isoparametric space-time FEM on
    moving domains. Master's thesis, NAM, University of Göttingen, 2018.
[3] R. I. Saye, High-order quadrature on multi-component domains implicitly
    defined by multivariate polynomials, J. Comput. Phys. 448:110720, 2022.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
import sys
from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh
from xfem import *

ngsglobals.msg_level = 1

# -------------------------------- PARAMETERS ---------------------------------
# Quadrilateral (or simplicial) mesh
quad_mesh = True
# Finite element space order (and order of the level set function)
order = 3
# Mesh sizes h = 2/n
n_levels = [8, 16, 32, 64]
# Stabilization parameter for ghost-penalty
gamma_stab = 0.1
# Stabilization parameter for Nitsche
lambda_nitsche = 10 * order * order

if len(sys.argv) > 1 and sys.argv[1] == "testmode":
    n_levels = [8, 16]

# Quadrature on cut elements: Algoim for higher order level sets (without mesh
# deformation), the standard rules for piecewise linear level sets
ngsxfemglobals.cut_quadrature = "auto"

# ----------------------------------- MAIN ------------------------------------
# Manufactured exact solution for monitoring the error
r2 = 3 / 4  # outer radius
r1 = 1 / 4  # inner radius
r = sqrt(x**2 + y**2)
levelset = (x**2 + y**2 - r1**2) * (x**2 + y**2 - r2**2)

exact = (20 * (r2 - r) * (r - r1)).Compile()
coeff_f = - (exact.Diff(x).Diff(x) + exact.Diff(y).Diff(y)).Compile()


def solve(n, lset_order):
    mesh = MakeStructured2DMesh(quads=quad_mesh, nx=n, ny=n,
                                mapping=lambda x, y: (2 * x - 1, 2 * y - 1))

    # level set function: interpolation into the space of order lset_order
    lset = GridFunction(H1(mesh, order=lset_order))
    lset.Set(levelset)

    # Element, facet and dof marking w.r.t. the (discrete) level set function
    ci = CutInfo(mesh, lset)
    hasneg = ci.GetElementsOfType(HASNEG)
    hasif = ci.GetElementsOfType(IF)

    # facets used for stabilization:
    ba_facets = GetFacetsWithNeighborTypes(mesh, a=hasneg, b=hasif)

    Vhbase = H1(mesh, order=order, dirichlet=[], dgjumps=True)
    Vh = Restrict(Vhbase, hasneg)

    gfu = GridFunction(Vh)

    u, v = Vh.TrialFunction(), Vh.TestFunction()
    h = specialcf.mesh_size
    nv = Normalize(grad(lset))

    # integration domains:
    dx = dCut(lset, NEG, definedonelements=hasneg)
    ds = dCut(lset, IF, definedonelements=hasif)
    dw = dFacetPatch(definedonelements=ba_facets)

    a = BilinearForm(Vh, symmetric=False)
    # Diffusion term
    a += grad(u) * grad(v) * dx
    # Nitsche term
    a += -grad(u) * nv * v * ds
    a += -grad(v) * nv * u * ds
    a += (lambda_nitsche / h) * u * v * ds
    # Ghost penalty stabilization (near the boundary)
    a += gamma_stab / h**2 * (u - u.Other()) * (v - v.Other()) * dw

    # R.h.s. term:
    f = LinearForm(Vh)
    f += coeff_f * v * dx

    # Assemble system
    a.Assemble()
    f.Assemble()

    # Solve linear system
    gfu.vec.data = a.mat.Inverse(Vh.FreeDofs()) * f.vec

    # Measure the error
    l2error = sqrt(Integrate((gfu - exact)**2 * dx, mesh))
    return mesh, lset, gfu, l2error


errors = {1: [], order: []}
print("L2 errors with level set of order 1 (standard rules) and of order "
      "{0} (Algoim) on {1} meshes:".format(order, "quadrilateral" if quad_mesh else "triangular"))
print("{:>5} {:>12} {:>6} {:>12} {:>6}".format("n", "P1-lset", "eoc", "P{}-lset".format(order), "eoc"))
for i, n in enumerate(n_levels):
    line = "{:5d}".format(n)
    for lset_order in errors:
        mesh, lset, gfu, err = solve(n, lset_order)
        errors[lset_order].append(err)
        eoc = log(errors[lset_order][-2] / err) / log(2) if i > 0 else 0
        line += " {:12.4e} {:6.2f}".format(err, eoc)
    print(line)

if len(sys.argv) > 1 and sys.argv[1] == "testmode":
    # higher order convergence with Algoim (and much more accurate than with P1)
    assert errors[order][-1] < 1e-3
    assert errors[order][-1] < 0.05 * errors[1][-1]
    assert errors[order][0] / errors[order][-1] > 2**(order + 0.5)

# visualization (finest mesh, higher order level set):
DrawDC(lset, gfu, 0, mesh, "uh")
