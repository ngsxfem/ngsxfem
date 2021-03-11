"""
In this example we integrate over a circle with radius 1 to illustrate
the higher order convergence achieved through the isoparametric mapping
approach. The mesh is deformed in order to map the P1 level set
approximation onto a higher-order approximation of the smooth level set,
thereby overcoming the geometry error of h^2 made by the piecewise
linear level set approximation and recovering the full integration
order.

Literature
----------
[1] C. Lehrenfeld. High order unfitted finite element methods on level
    set domains using isoparametric mappings. Comp. Meth. Appl. Mech.
    Eng., 300:716-733, 2016.
[2] C. Lehrenfeld. A higher order isoparametric fictitious domain method
    for level set domains. In S. P. A. Bordas, E. Burman, M. G. Larson,
    and M. A. Olshanskii, editors, Geometrically Unfitted Finite Element
    Methods and Applications, pages 65-92. Springer International
    Publishing, 2017.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *

from xfem import *
from xfem.lsetcurv import *

from math import pi


# -------------------------------- PARAMETERS ---------------------------------
maxh = 0.5
order = 5

maxreflvl = 6


# ----------------------------------- MAIN ------------------------------------
def Make2DProblem(maxh=2):
    square = SplineGeometry()
    square.AddRectangle([-1, -1], [1, 1], bc=1)
    mesh = Mesh(square.GenerateMesh(maxh=maxh, quad_dominated=False))
    return mesh


# Circle with radius 0.5
levelset = sqrt(x * x + y * y) - 0.5
integrand = CoefficientFunction(1.0)
referencevals = {POS: 4.0 - 0.25 * pi, NEG: 0.25 * pi, IF: pi}

mesh = Make2DProblem(maxh=maxh)

# Level set mesh adaptation machinery
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2,
                                      discontinuous_qn=True)
lsetp1 = lsetmeshadap.lset_p1

# Containers to compute error convergence
err_uncurved, err_curved = {}, {}
eoc_uncurved, eoc_curved = {}, {}

for key in [NEG, POS, IF]:
    err_curved[key] = []
    err_uncurved[key] = []
    eoc_curved[key] = []
    eoc_uncurved[key] = []


# Main loop
for reflevel in range(maxreflvl):
    if(reflevel > 0):
        mesh.Refine()
    lsetmeshadap.CalcDeformation(levelset)

    for key in [NEG, POS, IF]:
        dx = dCut(lsetp1, key, order=order)
        dy = dCut(lsetp1, key, order=order, deformation=lsetmeshadap.deform)
        integral_uncurved = Integrate(integrand*dx, mesh)
        integral_curved = Integrate(integrand*dy, mesh)

        err_curved[key].append(abs(integral_curved - referencevals[key]))
        err_uncurved[key].append(abs(integral_uncurved - referencevals[key]))

    # Mark cut elements for refinement:
    RefineAtLevelSet(gf=lsetmeshadap.lset_p1)

for key in [NEG, POS, IF]:
    eoc_curved[key] = [log(a / b) / log(2) for (a, b)
                       in zip(err_curved[key][:-1], err_curved[key][1:])]
    eoc_uncurved[key] = [log(a / b) / log(2) for (a, b)
                         in zip(err_uncurved[key][:-1], err_uncurved[key][1:])]

print("errors (uncurved):  \n{}\n".format(err_uncurved))
print("   eoc (uncurved):  \n{}\n".format(eoc_uncurved))
print("errors (  curved):  \n{}\n".format(err_curved))
print("   eoc (  curved):  \n{}\n".format(eoc_curved))

Draw(levelset, mesh, "levelset")
Draw(lsetmeshadap.deform, mesh, "deformation")
Draw(lsetmeshadap.lset_p1, mesh, "levelset(P1)")
