"""
    Test: Integrate on an unfitted unit square
"""
# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *

# -------------------------------- PARAMETERS ---------------------------------
h_max = 0.4
k = 4

function = x * (1 - x) * y * (1 - y)
value = 1 / 36

# ----------------------------------- DATA ------------------------------------


def level_sets():
    return [-y, x - 1, y - 1, -x]


nr_ls = len(level_sets())

# ------------------------------ BACKGROUND MESH ------------------------------
geo = SplineGeometry()
geo.AddRectangle((-0.2, -0.2), (1.2, 1.2),
                 bcs=("bottom", "right", "top", "left"))
mesh = Mesh(geo.GenerateMesh(maxh=h_max))

# ---------------------------- LEVELSET & CUT-INFO ----------------------------

level_sets_p1 = tuple(GridFunction(H1(mesh, order=1)) for i in range(nr_ls))
for i, lsetp1 in enumerate(level_sets_p1):
    InterpolateToP1(level_sets()[i], lsetp1)

square = DomainTypeArray([(NEG, NEG, NEG, NEG)])
Draw(square.Indicator(level_sets_p1), mesh, "mlsets")

lset_dom_inner = {"levelset": level_sets_p1, "domain_type": square}

# --------------------------------- INTEGRATE ---------------------------------
result = Integrate(levelset_domain=lset_dom_inner, mesh=mesh, cf=function, order=k)
print("Result: {}".format(result))
print("Error : {}".format(abs(result-value)))
assert abs(result-value) < 1e-12