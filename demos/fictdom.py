"""
In this example we solve a scalar *unfitted* PDE problem. As a
discretisation method we use a level set based geometry description and
a Cut (or Fictitious) Finite element method with a Nitsche formulation
to impose boundary conditions. For stability we add a ghost penalty
stabilization.

Domain:
-------
The domain is [-1,1]^2 while the interface is a ringe described by a
level set function. In the discretisation the level set function is
approximated with a piecewise linear interpolation. This approximate
geometry is then mapped by applying a mesh deformation resulting in
a higher order geometry approximation.

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
* Geometry approximation using isoparametric unfitted FEM

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
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.internal import *
from xfem import *
from xfem.lsetcurv import *

ngsglobals.msg_level = 2

# -------------------------------- PARAMETERS ---------------------------------
# Quadrilateral (or simplicial mesh)
quad_mesh = False
# Mesh diameter
maxh = 0.1
# Finite element space order
order = 3
# Stabilization parameter for ghost-penalty
gamma_stab = 0.1
# Stabilization parameter for Nitsche
lambda_nitsche = 10 * order * order


# ----------------------------------- MAIN ------------------------------------
# Geometry and Mesh
square = SplineGeometry()
square.AddRectangle((-1, -1), (1, 1), bc=1)
ngmesh = square.GenerateMesh(maxh=maxh, quad_dominated=quad_mesh)
mesh = Mesh(ngmesh)


# Manufactured exact solution for monitoring the error
r2 = 3 / 4  # outer radius
r1 = 1 / 4  # inner radius
rc = (r1 + r2) / 2.0
rr = (r2 - r1) / 2.0
r = sqrt(x**2 + y**2)
levelset = IfPos(r - rc, r - rc - rr, rc - r - rr)

exact = (20 * (r2 - sqrt(x**2 + y**2)) * (sqrt(x**2 + y**2) - r1)).Compile()
coeff_f = - (exact.Diff(x).Diff(x) + exact.Diff(y).Diff(y)).Compile()

# Higher order level set approximation
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.1,
                                      discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lsetp1 = lsetmeshadap.lset_p1


# Element, facet and dof marking w.r.t. boundary approximation with lsetp1:
ci = CutInfo(mesh, lsetp1)
hasneg = ci.GetElementsOfType(HASNEG)
hasif = ci.GetElementsOfType(IF)

# facets used for stabilization:
ba_facets = GetFacetsWithNeighborTypes(mesh, a=hasneg, b=hasif)

Vhbase = H1(mesh, order=order, dirichlet=[], dgjumps=True)
Vh = Restrict(Vhbase, hasneg)

gfu = GridFunction(Vh)

u, v = Vh.TrialFunction(), Vh.TestFunction()
h = specialcf.mesh_size
n = Normalize(grad(lsetp1))

# integration domains:
dx = dCut(lsetp1, NEG, definedonelements=hasneg, deformation=deformation)
ds = dCut(lsetp1, IF, definedonelements=hasif, deformation=deformation)
dw = dFacetPatch(definedonelements=ba_facets, deformation=deformation)

a = BilinearForm(Vh, symmetric=False)
# Diffusion term
a += grad(u) * grad(v) * dx
# Nitsche term
a += -grad(u) * n * v * ds
a += -grad(v) * n * u * ds
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
print("L2 Error: {0}".format(l2error))

# visualization:
Draw(deformation, mesh, "deformation")
DrawDC(lsetp1, gfu, 0, mesh, "uh", deformation=deformation)
