"""
In this example we integrate over a quarter circle on a series of
quadrilateral meshes.
"""


# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.lsetcurv import *
from math import pi


# -------------------------------- PARAMETERS ---------------------------------
# Radius of circle
r = 0.6
# Number of mesh refinements
n_ref = 8
# Integration order
order = 2


# ----------------------------------- MAIN ------------------------------------

square = SplineGeometry()
square.AddRectangle([0, 0], [1, 1], bc=1)
# mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=False))
mesh = Mesh(square.GenerateMesh(maxh=100, quad_dominated=True))

# domains = [NEG,POS,IF]
domains = [NEG, POS, IF]

levelset = sqrt(x * x + y * y) - r
referencevals = {POS: 1 - pi * r * r / 4, NEG: pi * r * r / 4, IF: r * pi / 2}


errors = dict()
eoc = dict()

for key in domains:
    errors[key] = []
    eoc[key] = []

for i in range(n_ref):
    lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2,
                                          discontinuous_qn=True)
    V = H1(mesh, order=1)
    lset_approx = GridFunction(V)
    InterpolateToP1(levelset, lset_approx)
    # Draw(lset_approx,mesh,"lset_p1")

    f = CoefficientFunction(1)

    deformation = lsetmeshadap.CalcDeformation(levelset)
    # mesh.SetDeformation(deformation)
    # Draw(deformation,mesh,"deformation")

    for key in domains:
        lset_dom = {"levelset": lset_approx, "domain_type": key}
        integral = Integrate(levelset_domain=lset_dom,
                             cf=f, mesh=mesh, order=order)
        print("Result of Integration Reflevel ",
              i, ", Key ", key, " : ", integral)
        errors[key].append(abs(integral - referencevals[key]))

    mesh.deformation = None

    if i < n_ref - 1:
        # RefineAtLevelSet(gf=lset_approx)
        mesh.Refine()

for key in domains:
    eoc[key] = [log(errors[key][i + 1] / errors[key][i]) / log(0.5)
                for i in range(n_ref - 1)]

print("L2-errors:", errors)
# l2_eoc = [log(errors[i+1]/errors[i])/log(0.5) for i in range(n_ref-1)]
print("experimental order of convergence (L2):", eoc)

for key in domains:
    for i in range(2, len(eoc[key])):
        eoc_v = eoc[key][i]
        if eoc_v < 1.8:
            raise RuntimeError("Order of Convergence < 1.8 detected!")
