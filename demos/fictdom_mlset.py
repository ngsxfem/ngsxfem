"""
In this example we solve an unfitted Poisson problem similar to the one in
`fictdom.py`, however this time with the unfitted geometry being the
unit square. This example shall illustrate the functionality of ngsxfem to
solve PDE problems on geometries described via multiple level set functions.

PDE problem + Discretisation + Geometry + Implementation aspects:
-----------------------------------------------------------------
* As in fictdom.py except for the different geometry and its handling.

Used Features:
--------------
* Quadrature with respect to multiple level set functions., see the
  'mlset_pde' jupyter tutorial.

* MultiLevelsetCutInfo, see the 'mlset_basic' jupyter tutorial.

* DomainTypeArray convenience layer, see the 'mlset_basic' jupyter
  tutorial.

* Restricted BilinearForm, jupyter tutorial `basics`.

* Cut Differential Symbols, jupyter tutorials `intlset` and `cutfem`.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *

ngsglobals.msg_level = 2

# -------------------------------- PARAMETERS ---------------------------------
# Domain corners
ll, ur = (-0.2, -0.2), (1.2, 1.2)
# Initial mesh diameter
initial_maxh = 0.4
# Number of mesh bisections
nref = 3
# Order of finite element space
k = 1

# Stabilization parameter for ghost-penalty
gamma_s = 0.5
# Stabilization parameter for Nitsche
gamma_n = 10

# ----------------------------------- MAIN ------------------------------------
# Set up the level sets, exact solution and right-hand side


def level_sets():
    return [-y, x - 1, y - 1, -x]


nr_ls = len(level_sets())
u_ex = 16 * x * (1 - x) * y * (1 - y)
grad_u_ex = (u_ex.Diff(x).Compile(), u_ex.Diff(y).Compile())
rhs = -(u_ex.Diff(x).Diff(x) + u_ex.Diff(y).Diff(y)).Compile()

# Geometry and mesh
geo = SplineGeometry()
geo.AddRectangle(ll, ur, bcs=("bottom", "right", "top", "left"))
ngmesh = geo.GenerateMesh(maxh=initial_maxh)
for i in range(nref):
    ngmesh.Refine()
mesh = Mesh(ngmesh)


# Level set and cut-information
P1 = H1(mesh, order=1)
lsetsp1 = tuple(GridFunction(P1) for i in range(nr_ls))
for i, lsetp1 in enumerate(lsetsp1):
    InterpolateToP1(level_sets()[i], lsetp1)
    Draw(lsetp1, mesh, "lsetp1_{}".format(i))

square = DomainTypeArray((NEG, NEG, NEG, NEG))
with TaskManager():
    square.Compress(lsetsp1)
    boundary = square.Boundary()
    boundary.Compress(lsetsp1)

mlci = MultiLevelsetCutInfo(mesh, lsetsp1)


# Element and degrees-of-freedom markers
els_if_singe = {dtt: BitArray(mesh.ne) for dtt in boundary}
facets_gp = BitArray(mesh.nedge)

hasneg = mlci.GetElementsWithContribution(square)

# Finite element space
Vhbase = H1(mesh, order=k, dgjumps=True)
Vh = Restrict(Vhbase, hasneg)
gfu = GridFunction(Vh)

hasif = mlci.GetElementsWithContribution(boundary)
Draw(BitArrayCF(hasif), mesh, "hasif")

for i, (dtt, els_bnd) in enumerate(els_if_singe.items()):
    els_bnd[:] = mlci.GetElementsWithContribution(dtt)
    Draw(BitArrayCF(els_bnd), mesh, "els_if_singe" + str(i))

facets_gp = GetFacetsWithNeighborTypes(mesh, a=hasneg, b=hasif,
                                       use_and=True)

els_gp = GetElementsWithNeighborFacets(mesh, facets_gp)
Draw(BitArrayCF(els_gp), mesh, "gp_elements")

# Bilinear and linear forms of the weak formulation
u, v = Vh.TnT()
h = specialcf.mesh_size
normals = square.GetOuterNormals(lsetsp1)

# Set up the integrator symbols
dx = dCut(lsetsp1, square, definedonelements=hasneg)
ds = {dtt: dCut(lsetsp1, dtt, definedonelements=els_if_singe[dtt])
      for dtt in boundary}
dw = dFacetPatch(definedonelements=facets_gp)

# Construct integrator
a = RestrictedBilinearForm(Vh, facet_restriction=facets_gp, check_unused=False)
a += InnerProduct(grad(u), grad(v)) * dx
for bnd, n in normals.items():
    a += -InnerProduct(grad(u) * n, v) * ds[bnd]
    a += -InnerProduct(grad(v) * n, u) * ds[bnd]
    a += (gamma_n * k * k / h) * InnerProduct(u, v) * ds[bnd]
a += gamma_s / (h**2) * (u - u.Other()) * (v - v.Other()) * dw

f = LinearForm(Vh)
f += rhs * v * dx


# Assemble and solve the linear system
f.Assemble()
a.Assemble()

gfu.vec.data = a.mat.Inverse(Vh.FreeDofs()) * f.vec

Draw(gfu, mesh, "uh")

# Post-processing
err_l2 = sqrt(Integrate((gfu - u_ex)**2 * dx.order(2 * k), mesh))
err_h1 = sqrt(Integrate((Grad(gfu) - grad_u_ex)**2 * dx.order(2 * (k - 1)),
                        mesh))

print("L2 error = {:1.5e}".format(err_l2), "H1 error = {:1.5e}".format(err_h1))
