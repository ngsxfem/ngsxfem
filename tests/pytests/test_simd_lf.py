from netgen.geom2d import SplineGeometry
from xfem import ngsxfemglobals
import numpy as np
from math import pi
from timeit import default_timer as timer
import pytest
from ngsolve import *
from xfem import *
from xfem.mlset import *


def level_sets():
    return [-y, x - 1, y - 1, -x]


rhs = 32 * (y * (1 - y) + x * (1 - x))


h0 = 0.05
k = 3
levelsets = level_sets()
nr_ls = len(levelsets)
geo = SplineGeometry()
geo.AddRectangle((-0.2, -0.2), (1.2, 1.2),
                 bcs=("bottom", "right", "top", "left"))

mesh = Mesh(geo.GenerateMesh(maxh=0.5))
# ------------------------- Finite Element Space --------------------------
V = H1(mesh, order=k, dgjumps=True)

gfu = GridFunction(V)
freedofs = BitArray(V.ndof)

# -------------------------- Levelset & Cut-Info --------------------------
level_sets_p1 = tuple(GridFunction(H1(mesh, order=1))
                      for i in range(nr_ls))
for i, lsetp1 in enumerate(level_sets_p1):
    InterpolateToP1(levelsets[i], lsetp1)

square = DomainTypeArray([(NEG, NEG, NEG, NEG)])
square.Compress(level_sets_p1)
lset_dom_inner = {"levelset": level_sets_p1, "domain_type": square}

boundary = square.Boundary()
lsets_bnd = {}
for dtt in boundary:
    lsets_bnd[dtt] = {"levelset": level_sets_p1, "domain_type": dtt}


# --------------------------- (Bi)Linear Forms ----------------------------
u, v = V.TnT()


forcing = rhs * v

ngsxfemglobals.SetDefaults()

f = LinearForm(V)
f += SymbolicLFI(lset_dom_inner, form=forcing)

g = LinearForm(V)
g += SymbolicLFI(lset_dom_inner, form=forcing)


start = timer()

f.Assemble()
end = timer()
vals1 = f.vec.FV().NumPy()

print("Elapsed time: ")
print(end-start)
ngsxfemglobals.SwitchSIMD(True)
start = timer()

g.Assemble()
end = timer()
vals2 = g.vec.FV().NumPy()
print("Elapsed time: ")
print(end-start)
#assert(np.linalg.norm(vals1-vals2) < 1e-8)
print(np.linalg.norm(vals1-vals2))
print(vals2-vals1)
# print(f.vec.data)
# ----------------------------- Solve Problem -----------------------------
