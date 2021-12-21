from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem import ngsxfemglobals
from xfem.lsetcurv import LevelSetMeshAdaptation
import numpy as np
from math import pi
from timeit import default_timer as timer

# -------------------------------- PARAMETERS ---------------------------------
# Domain corners
ll, ur = (-1.5, -1.5), (1.5, 1.5)
# Mesh size
maxh = 0.2
# Finite element space order
order = 2

# Diffusion coefficients for the sub-domains (NEG/POS):
alpha = [1.0, 2.0]
# Nitsche penalty parameter
lambda_nitsche = 20

formulation = "XFEM"  # or "CUTFEM"

# ----------------------------------- MAIN ------------------------------------
# We generate the background mesh of the domain and use a simplicity
# triangulation to obtain a mesh with quadrilaterals use
# 'quad_dominated=True'
square = SplineGeometry()
square.AddRectangle(ll, ur, bc=1)
mesh = Mesh(square.GenerateMesh(maxh=maxh, quad_dominated=False))

r44 = x**4 + y**4

r41 = sqrt(sqrt(r44))
solution = [1 + pi / 2 - sqrt(2.0) * cos(pi / 4 * r44), pi / 2 * r41]
coef_f = [-alpha[i] * (solution[i].Diff(x).Diff(x)
                    + solution[i].Diff(y).Diff(y)) for i in range(2)]

# Level set function of the domain (phi = ||x||_4 - 1) and its interpolation:
levelset = r41 - 1.0


lsetadap = NoDeformation(mesh, levelset)
lsetp1 = lsetadap.lset_p1

# Background FESpaces (used as CutFESpaces later-on):
Vh = H1(mesh, order=order, dirichlet=".*")

# Gathering information on cut elements:
ci = CutInfo(mesh, lsetp1)

Vhx = XFESpace(Vh, lsetp1)
VhG = Vh * Vhx


(u_std, u_x), (v_std, v_x) = VhG.TnT()

u = [u_std + op(u_x) for op in [neg, pos]]
v = [v_std + op(v_x) for op in [neg, pos]]


# Integration domains for integration on negative/positive sub-domains
# and on the interface: Here, the integration is (geometrically) exact
# if the "levelset"-argument is a piecewise (multi-)linear function.
# We further provide a mesh deformation that is applied in the higher order
# case:
ngsxfemglobals.SetDefaults()

dx = tuple([dCut(lsetp1, dt, deformation=lsetadap.deform,
                definedonelements=ci.GetElementsOfType(HAS(dt)))
            for dt in [NEG, POS]])

# R.h.s.:
f = LinearForm(VhG)
f += sum(coef_f[i] * v[i] * dx[i] for i in [0, 1])
start =timer()

f.Assemble()
end = timer()
vals1 = f.vec.FV().NumPy()
print(vals1)
print("Elapsed time: ")
print(end-start)
ngsxfemglobals.simd_eval = True
start =timer()

f.Assemble()
end = timer()
vals2 = f.vec.FV().NumPy()
print("Elapsed time: ")
print(end-start)
assert(np.linalg.norm(vals1-vals2)<1e-8)
#print(f.vec.data)