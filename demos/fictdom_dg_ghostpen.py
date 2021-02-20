"""
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.internal import *
from xfem import *
from xfem.lsetcurv import *


# -------------------------------- PARAMETERS ---------------------------------
quad_mesh = False
maxh = 0.1
order = 2

lambda_nitsche = 10 * order * order
lambda_dg = 10 * order * order


# ----------------------------------- MAIN ------------------------------------
# Geometry and Mesh
square = SplineGeometry()
square.AddRectangle((-1, -1), (1, 1), bc=1)

ngmesh = square.GenerateMesh(maxh=maxh, quad_dominated=quad_mesh)
mesh = Mesh(ngmesh)


# Manufactured exact solution
r2 = 3 / 4  # outer radius
r1 = 1 / 4  # inner radius
rc = (r1 + r2) / 2.0
rr = (r2 - r1) / 2.0
r = sqrt(x**2 + y**2)
levelset = IfPos(r - rc, r - rc - rr, rc - r - rr)

coeff_f = CoefficientFunction(-20 * ((r1 + r2) / sqrt(x**2 + y**2) - 4))

# for monitoring the error
exact = CoefficientFunction(
    20 * (r2 - sqrt(x**2 + y**2)) * (sqrt(x**2 + y**2) - r1)).Compile()


# Higher order level set approximation
lsetmeshadap = LevelSetMeshAdaptation(
    mesh, order=order, threshold=0.1, discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lsetp1 = lsetmeshadap.lset_p1

lset_neg = {"levelset": lsetp1, "domain_type": NEG, "subdivlvl": 0}
lset_if = {"levelset": lsetp1, "domain_type": IF, "subdivlvl": 0}


# Element, facet and dof marking w.r.t. boundary approximation with lsetp1:
ci = CutInfo(mesh, lsetp1)
hasneg = ci.GetElementsOfType(HASNEG)

Vh = L2(mesh, order=order, dirichlet=[], dgjumps=True)
active_dofs = GetDofsOfElements(Vh, hasneg)
Vh = Compress(Vh, active_dofs)

gfu = GridFunction(Vh)

hasif = ci.GetElementsOfType(IF)

ba_gp_facets = GetFacetsWithNeighborTypes(mesh, a=hasneg, b=hasif)
ba_fd_facets = GetFacetsWithNeighborTypes(mesh, a=hasneg, b=hasneg)

a = RestrictedBilinearForm(Vh, "a", hasneg, ba_fd_facets, check_unused=False)
f = LinearForm(Vh)

u, v = Vh.TnT()
h = specialcf.mesh_size
nh = 1.0 / Norm(grad(lsetp1)) * grad(lsetp1)

# element-wise diffusion term
dx = dCut(lsetp1, NEG, definedonelements=hasneg, deformation=deformation)
a += grad(u) * grad(v) * dx

# ghost penalty
dw = dFacetPatch(definedonelements=ba_gp_facets, deformation=deformation)
a += 0.1 / h**2 * (u - u.Other()) * (v - v.Other()) * dw

# Interior penalty terms:
nF = specialcf.normal(mesh.dim)
flux_u = -0.5 * (grad(u) + grad(u.Other())) * nF
flux_v = -0.5 * (grad(v) + grad(v.Other())) * nF
jump_u = u - u.Other()
jump_v = v - v.Other()
dk = dCut(lsetp1, NEG, skeleton=True, definedonelements=ba_fd_facets,
          deformation=deformation)
a += (lambda_dg / h * jump_u * jump_v + flux_u * jump_v + flux_v * jump_u) * dk

# Nitsche term
ds = dCut(lsetp1, IF, definedonelements=hasif, deformation=deformation)
a += (-grad(u)*nh * v - grad(v)*nh * u + lambda_nitsche/h * u * v) * ds

# R.h.s. term:
f += coeff_f * v * dx

# Apply mesh adaptation

# Assemble system
f.Assemble()
mesh.deformation = deformation
a.Assemble()

# Solve linear system
gfu.vec.data = a.mat.Inverse(Vh.FreeDofs(), "sparsecholesky") * f.vec

# measure the error
l2error = sqrt(Integrate(lset_neg, (gfu - exact) * (gfu - exact), mesh))
print("L2 Error: {0}".format(l2error))

# Unset mesh adaptation
mesh.deformation = None


# Visualization:
Draw(deformation, mesh, "deformation")
Draw(levelset, mesh, "levelset")
Draw(lsetp1, mesh, "lsetp1")
Draw(gfu, mesh, "extu")
Draw(IfPos(-lsetp1, gfu, float('nan')), mesh, "u")
warped_u = CoefficientFunction((deformation[0],
                                deformation[1],
                                IfPos(-lsetp1, 0.2 * gfu, float('nan'))))
Draw(warped_u, mesh, "warped_u")

visoptions.autoscale = False
visoptions.mminval = 0
visoptions.mmaxval = 1.25
visoptions.deformation = 1
