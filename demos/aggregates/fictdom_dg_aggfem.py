"""
In this example we solve a scalar *unfitted* PDE problem as in fictdom.py
but with an unfitted DG discretization.

Domain + PDE problem:
---------------------
As in fictdom.py

Discretisation:
---------------
* Background discontinuous FEspace restricted to active domain (CutFEM)

* Interior penalty DG formulation inside the domain, see e.g. [1]

* Nitsche formulation to impose boundary conditions, see. e.g. [2]

* Element Aggregation to deal with bad cuts.

Implementational aspects:
-------------------------
* Geometry approximation using isoparametric unfitted FEM

* A (sparse) direct solver is applied to solve the arising linear systems.

References:
-----------
* All concepts that are used here are explained in the jupyter tutorials
  `basics`, `intlset` and `cutfem`.


Literature:
-----------
[1] C. Gürkan, S. Sticko, A. Massing, Stabilized Cut Discontinuous Galerkin
    Methods for Advection-Reaction Problems, SIAM J. Sci. Comp. 42(3):
    A2620-A2654, 2020.
[2] E. Burman, P. Hansbo, Fictitious domain finite element methods using
    cut elements: II. A stabilized Nitsche method, Appl. Num. Math.
    62(4):328-341, 2012.
[3] S. Badia, F. Verdugo, and A. F. Martín. The aggregated unfitted finite 
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
maxh = 0.1
# Finite element space order
order = 2

# Stabilization parameter for Nitsche
lambda_nitsche = 10 * order * order
# Stabilization parameter for interior penalty
lambda_dg = 10 * order * order

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
# Facets with parts inside the domain
ba_fd_facets = GetFacetsWithNeighborTypes(mesh, a=els_hasneg, b=els_hasneg)

Vhbase = L2(mesh, order=order, dirichlet=[], dgjumps=True)
Vh = Restrict(Vhbase, els_hasneg)

gfu = GridFunction(Vh)

u, v = Vh.TnT()
h = specialcf.mesh_size
nh = 1.0 / Norm(grad(lsetp1)) * grad(lsetp1)
nF = specialcf.normal(mesh.dim)
flux_u = -0.5 * (grad(u) + grad(u.Other())) * nF
flux_v = -0.5 * (grad(v) + grad(v.Other())) * nF
jump_u = u - u.Other()
jump_v = v - v.Other()


# Element-wise integrals
dx = dCut(lsetp1, NEG, definedonelements=els_hasneg, deformation=deformation)
# Interior skeleton integrals:
dk = dCut(lsetp1, NEG, skeleton=True, definedonelements=ba_fd_facets,
          deformation=deformation)
# Domain boundary integrals
ds = dCut(lsetp1, IF, definedonelements=els_if, deformation=deformation)

a = RestrictedBilinearForm(Vh, "a", els_hasneg, ba_fd_facets, check_unused=False)
a += grad(u) * grad(v) * dx
a += (lambda_dg / h * jump_u * jump_v + flux_u * jump_v + flux_v * jump_u) * dk
a += (-grad(u) * nh * v - grad(v) * nh * u + lambda_nitsche / h * u * v) * ds

f = LinearForm(Vh)
f += coeff_f * v * dx

# Assemble system
f.Assemble()
a.Assemble()

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

# Unset mesh adaptation
mesh.deformation = None

# Visualization:
Draw(deformation, mesh, "deformation")
DrawDC(lsetp1, gfu, 0, mesh, "uh")
