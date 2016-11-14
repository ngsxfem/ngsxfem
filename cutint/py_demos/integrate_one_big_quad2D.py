# integration on lset domains

from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *

from netgen.geom2d import SplineGeometry

square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
#mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=False))
mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))

levelset = 1-2*x-2*y
referencevals = { POS : 1./8, NEG : 7./8, IF : 1./sqrt(2)}

V = H1(mesh,order=1)
lset_approx = GridFunction(V)
InterpolateToP1(levelset,lset_approx)

f = CoefficientFunction (1.0)

domains = [NEG,POS,IF]

errors = dict()
eoc = dict()

for key in domains:
    errors[key] = []
    eoc[key] = []

for order in range(10):
    for key in domains:
        integral = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=f,order=order,domain_type=key,heapsize=1000000)
        print("Integral on Domain ", key, " : ",integral)
        errors[key].append(abs(integral - referencevals[key]))

Draw(levelset, mesh, "lset")
print("L2 Errors:", errors)
