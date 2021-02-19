"""
In this example we compute a Poisson problem on the unit square to
illustrate the functionality of ngsxfem to solve PDE problems on
geometries described via multiple level sets.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.solvers import PreconditionedRichardson as PreRic
from xfem import *
from xfem.mlset import *

ngsglobals.msg_level = 2
SetNumThreads(4)


# -------------------------------- PARAMETERS ---------------------------------
ll, ur = (-0.2, -0.2), (1.2, 1.2)
h0 = 0.4
Lx = 3
k = 1

gamma_n = 10
gamma_s = 0.5

inverse = "sparsecholesky"


# ----------------------------------- MAIN ------------------------------------
# Set up the level sets, exact solution and right-hand side
def level_sets():
    return [-y, x - 1, y - 1, -x]


nr_ls = len(level_sets())
u_ex = 16 * x * (1 - x) * y * (1 - y)
grad_u_ex = CoefficientFunction((16 * (1 - 2 * x) * y * (1 - y),
                                 16 * x * (1 - x) * (1 - 2 * y)))
rhs = 32 * (y * (1 - y) + x * (1 - x))


# Geometry and mesh
geo = SplineGeometry()
geo.AddRectangle(ll, ur, bcs=("bottom", "right", "top", "left"))
ngmesh = geo.GenerateMesh(maxh=h0)
for i in range(Lx):
    ngmesh.Refine()
mesh = Mesh(ngmesh)

# Finite element space
V = H1(mesh, order=k, dgjumps=True)

gfu = GridFunction(V)
freedofs = BitArray(V.ndof)


# Level set and cut-information
P1 = H1(mesh, order=1)
level_sets_p1 = tuple(GridFunction(P1) for i in range(nr_ls))
for i, lsetp1 in enumerate(level_sets_p1):
    InterpolateToP1(level_sets()[i], lsetp1)
    Draw(lsetp1, mesh, "lsetp1_{}".format(i))

square = DomainTypeArray((NEG, NEG, NEG, NEG))
with TaskManager():
    square.Compress(level_sets_p1)
    boundary = square.Boundary()
    boundary.Compress(level_sets_p1)

mlci = MultiLevelsetCutInfo(mesh, level_sets_p1)


# Element and degrees-of-freedom markers
els_hasneg, els_if = BitArray(mesh.ne), BitArray(mesh.ne)
els_if_singe = {dtt: BitArray(mesh.ne) for dtt in boundary}
facets_gp = BitArray(mesh.nedge)

els_hasneg[:] = False
els_hasneg |= mlci.GetElementsWithContribution(square)

els_if[:] = False
els_if |= mlci.GetElementsWithContribution(boundary)
Draw(BitArrayCF(els_if), mesh, "els_if")

for i, (dtt, els_bnd) in enumerate(els_if_singe.items()):
    els_bnd[:] = False
    els_bnd |= mlci.GetElementsWithContribution(dtt)
    Draw(BitArrayCF(els_bnd), mesh, "els_if_singe" + str(i))

facets_gp[:] = False
facets_gp |= GetFacetsWithNeighborTypes(mesh, a=els_hasneg, b=els_if,
                                        use_and=True)

els_gp = GetElementsWithNeighborFacets(mesh, facets_gp)
Draw(BitArrayCF(els_gp), mesh, "gp_elements")

freedofs[:] = False
freedofs |= GetDofsOfElements(V, els_hasneg) & V.FreeDofs()


# Bilinear and linear forms of the weak formulation
u, v = V.TnT()
h = specialcf.mesh_size
normals = square.GetOuterNormals(level_sets_p1)

# Set up the integrator symbols
dx = dCut(level_sets_p1, square, definedonelements=els_hasneg)
ds = {dtt: dCut(level_sets_p1, dtt, definedonelements=els_if_singe[dtt])
      for dtt in boundary}
dw = dFacetPatch(definedonelements=facets_gp, skeleton=False)

# Construct integrator
a = RestrictedBilinearForm(V, element_restriction=els_hasneg,
                           facet_restriction=facets_gp, check_unused=False)
a += InnerProduct(Grad(u), Grad(v)) * dx
for bnd, n in normals.items():
    a += -InnerProduct(Grad(u) * n, v) * ds[bnd]
    a += -InnerProduct(Grad(v) * n, u) * ds[bnd]
    a += (gamma_n * k * k / h) * InnerProduct(u, v) * ds[bnd]
a += gamma_s / (h**2) * (u - u.Other()) * (v - v.Other()) * dw

f = LinearForm(V)
f += rhs * v * dx


# Assemble and solve the linear system
with TaskManager():
    f.Assemble()
    a.Assemble()
    inv = a.mat.Inverse(freedofs=freedofs, inverse=inverse)

    gfu.vec.data = PreRic(a=a, rhs=f.vec, pre=inv, freedofs=freedofs)

Draw(gfu, mesh, "solution")
Draw((gfu - u_ex), mesh, "err")


# Post-processing
with TaskManager():
    dx = dCut(level_sets_p1, square, order=2 * k)
    err_l2 = sqrt(Integrate((gfu - u_ex)**2 * dx, mesh))
    err_h1 = sqrt(Integrate((Grad(gfu) - grad_u_ex)**2 * dx, mesh))

print("\n")
print("L2 error = {:1.3e}".format(err_l2))
print("H1 error = {:1.3e}".format(err_h1))
