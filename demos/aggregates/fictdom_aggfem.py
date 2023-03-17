"""
In this example we solve a scalar *unfitted* PDE problem. As a
discretisation method we use a level set based geometry description and
a Cut (or Fictitious) Finite element method with a Nitsche formulation
to impose boundary conditions. For stability we use element aggregations
which constrain dofs on cut elements by extrapolating interior dofs

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

* Element aggregation to deal with bad cuts.

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
[2] S. Badia, F. Verdugo, and A. F. Martín. The aggregated unfitted finite 
    element method for elliptic problems. Comput. Methods Appl. Mech. 
    Engrg., 336:533–553, 2018.  
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
maxh = 0.1 / 2
# Finite element space order
order = 2
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
els_hasneg = ci.GetElementsOfType(HASNEG)
els_neg = ci.GetElementsOfType(NEG)
els_if = ci.GetElementsOfType(IF)

# Set up FE-Space
Vhbase = H1(mesh, order=order, dirichlet=[], dgjumps=False)
Vh = Restrict(Vhbase, els_hasneg)

gfu = GridFunction(Vh)

u, v = Vh.TrialFunction(), Vh.TestFunction()
h = specialcf.mesh_size
n = Normalize(grad(lsetp1))

# integration domains:
dx = dCut(lsetp1, NEG, definedonelements=els_hasneg, deformation=deformation)
ds = dCut(lsetp1, IF, definedonelements=els_if, deformation=deformation)

a = BilinearForm(Vh, symmetric=True)
# Diffusion term
a += grad(u) * grad(v) * dx
# Nitsche term
a += -grad(u) * n * v * ds
a += -grad(v) * n * u * ds
a += (lambda_nitsche / h) * u * v * ds

# R.h.s. term:
f = LinearForm(Vh)
f += coeff_f * v * dx

# Assemble unstable system
a.Assemble()
f.Assemble()

# Compute element patches for element aggregation
EA = ElementAggregation(mesh, els_neg, els_if)

# Set up embedding matrix for element aggregation
P = AggEmbedding(EA, Vh, deformation=deformation)

print(f'P shape: {P.shape}')
PT = P.CreateTranspose()

Aemb = PT @ a.mat @ P

# Solve linear system
gfuE = Aemb.Inverse(inverse='sparsecholesky') * (PT * f.vec)
gfu.vec.data = P * gfuE

# Measure the error
l2error = sqrt(Integrate((gfu - exact)**2 * dx, mesh))
print("L2 Error: {0}".format(l2error))

# visualization:
Draw(deformation, mesh, "deformation")
DrawDC(lsetp1, gfu, 0, mesh, "uh", deformation=deformation)
