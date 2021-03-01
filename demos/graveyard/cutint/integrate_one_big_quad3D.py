"""
integration on levelset domains
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve import *
from xfem import *


# ----------------------------------- MAIN ------------------------------------
cube = OrthoBrick(Pnt(0, 0, 0), Pnt(1, 1, 1)).bc(1)
geom = CSGeometry()
geom.Add(cube)
mesh = Mesh(geom.GenerateMesh(maxh=1, quad_dominated=True))


def binary_pow(x, a):
    if a == 0:
        return 0 * x + 1
    elif a == 1:
        return x
    else:
        print("Invalid argument a")


# levelset coefficients
c = [[[1, -2], [-2, 0]], [[-2, 0], [0, 0]]]

# levelset = 1-2*x-2*y-2*z #(sqrt(x*x+y*y+z*z)-0.5)
levelset = sum([c[alpha][beta][gamma] * binary_pow(x, alpha)
                * binary_pow(y, beta) * binary_pow(z, gamma)
                for alpha in [0, 1] for beta in [0, 1] for gamma in [0, 1]])
referencevals = {POS: 1. / 48, NEG: 47. / 48, IF: sqrt(3) / 8}

V = H1(mesh, order=1)
lset_approx = GridFunction(V)
InterpolateToP1(levelset, lset_approx)

f = CoefficientFunction(1.0)

domains = [NEG, POS, IF]

errors = dict()
eoc = dict()

for key in domains:
    errors[key] = []
    eoc[key] = []

for order in range(8):
    for key in domains:
        lset_dom = {"levelset": lset_approx, "domain_type": key}
        integral = Integrate(levelset_domain=lset_dom,
                             cf=f, mesh=mesh, order=order)
        print("Integral on Domain ", key, " : ", integral)
        errors[key].append(abs(integral - referencevals[key]))

Draw(levelset, mesh, "lset")
print("L2 Errors:", errors)
