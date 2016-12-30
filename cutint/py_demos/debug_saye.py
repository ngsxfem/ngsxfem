# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *

print("PR: \t",DebugSaye(-1, True))
print("PR: \t",DebugSaye(+1, True))
print("PR: \t",DebugSaye(-1, False))
print("PR: \t",DebugSaye(+1, False))

exit() #referencevals:  {DOMAIN_TYPE.POS: 0.817400277201530, DOMAIN_TYPE.NEG: 0.182599722798470}
print("\n\n----Second method-----\n\n")

from netgen.geom2d import SplineGeometry
square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))

lsetvals_list = [[1,-1,-3,-1]]

f = lambda x,y: 1.+0*x+0*y
f_ngs = f(x,y)
V = H1(mesh,order=1)
lset_approx = GridFunction(V)

domains = [NEG,POS]
levelset = 1-2*x-2*y
InterpolateToP1(levelset,lset_approx)

for key in domains:
    integral = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=f_ngs,order=1,domain_type=key,heapsize=1000000,use_saye = True)
    print(key," : ",integral)
