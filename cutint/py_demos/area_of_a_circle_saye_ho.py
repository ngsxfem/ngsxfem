from math import pi
from ngsolve import *
from xfem import *
from netgen.geom2d import SplineGeometry

square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
#mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=False))
mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))
r=0.6

domains = [NEG, POS, IF]

order = 3

levelset = sqrt(x*x+y*y)-r
referencevals = { POS : 1-pi*r*r/4, NEG : pi*r*r/4, IF : r*pi/2}

V = L2(mesh,order=order)
lset_approx = GridFunction(V)
#InterpolateToP1(levelset,lset_approx)
lset_approx.Set(levelset)
Draw(lset_approx,mesh,"lset_approx")

f = CoefficientFunction(1)

for key in domains:
    integral = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=f,order=order,domain_type=key,heapsize=1000000, use_saye = True)
    print("Result of Integration Key ",key," : ", integral)
    print("\t\tError: ", abs(integral - referencevals[key]))
