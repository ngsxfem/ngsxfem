# integration on lset domains

from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *

from netgen.csg import CSGeometry, OrthoBrick, Pnt
cube = OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) ).bc(1)
geom = CSGeometry()
geom.Add (cube)
ngmesh = geom.GenerateMesh(maxh=1, quad_dominated=True)
mesh = Mesh(ngmesh)

levelset = 1-2*x-2*y-2*z #(sqrt(x*x+y*y+z*z)-0.5)
referencevals = { POS : 1./48, NEG : 47./48, IF : sqrt(3)/8 }

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

for order in range(8):
    for key in domains:
        integral = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=f,order=order,domain_type=key,heapsize=1000000)
        print("Integral on Domain ", key, " : ",integral)
        errors[key].append(abs(integral - referencevals[key]))

Draw(levelset, mesh, "lset")
print("L2 Errors:", errors)
