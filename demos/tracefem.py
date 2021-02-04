"""
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.internal import *
from xfem import *


# -------------------------------- PARAMETERS ---------------------------------
# MEsh diameter
maxh = 0.3
# Polynomial order of FE space
order = 1
# Extra subdivisions of cut elements to generate quadrature rule
subdivlvl = 0


# ----------------------------------- MAIN ------------------------------------
# geometry
square = SplineGeometry()
square.AddRectangle([-1.5, -1.5], [1.5, 1.5], bc=1)
mesh = Mesh(square.GenerateMesh(maxh=maxh, quad_dominated=False))

levelset = sqrt(x * x + y * y) - 0.7

lset_approx = GridFunction(H1(mesh, order=1))
InterpolateToP1(levelset, lset_approx)

# Extended FESpace
VhG = H1(mesh, order=order, dirichlet=[])

# Overwrite freedofs of VhG to mark only dofs that are involved in the
# cut problem
ci = CutInfo(mesh, lset_approx)
ba_IF = ci.GetElementsOfType(IF)
cf_IF = BitArrayCF(ba_IF)
freedofs = VhG.FreeDofs()
freedofs &= GetDofsOfElements(VhG, ba_IF)

gfu = GridFunction(VhG)

# Coefficients / parameters:
n = 1.0 / Norm(grad(lset_approx)) * grad(lset_approx)
h = specialcf.mesh_size


# Tangential projection
def P(u):
    return u - (u * n) * n


# Expressions of test and trial functions:
u, v = VhG.TnT()

# Integration domains (and integration parameter "subdivlvl" and
# "force_intorder")
lset_if = {"levelset": lset_approx, "domain_type": IF, "subdivlvl": subdivlvl}

# Bilinear forms:
a = BilinearForm(VhG, symmetric=True)
a += SymbolicBFI(levelset_domain=lset_if, form=P(grad(u)) * P(grad(v)) + u * v)
a += SymbolicBFI(form=1.0 / h * (InnerProduct(grad(u), n)
                                 * InnerProduct(grad(v), n)),
                 definedonelements=ba_IF)
a.Assemble()

f = LinearForm(VhG)
f += SymbolicLFI(levelset_domain=lset_if, form=sin(x) * v)
f.Assemble()

gfu.vec[:] = 0.0
gfu.vec.data = a.mat.Inverse(freedofs) * f.vec

nan = CoefficientFunction(float('nan'))
Draw(IfPos(cf_IF - 0.5, gfu, nan), mesh, "u")

visoptions.mminval = -0.2
visoptions.mmaxval = 0.2
visoptions.deformation = 1
visoptions.autoscale = 0
