from math import pi
from ngsolve import *
from xfem.basics import *
from netgen.geom2d import SplineGeometry

square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
#mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=False))
mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))

#levelset values on the 4 vertices of the quad
#lsetvals = [1,0.9999,-1,-1] #case b)
#lsetvals = [1,2/3,-1,-2/3] #case c)
#lsetvals = [1,-1,1/3,-1] #case d)
#levelset = lsetvals[0]+(lsetvals[1] - lsetvals[0])*x +(lsetvals[3] - lsetvals[0])*y + (lsetvals[2]-lsetvals[1]-lsetvals[3]+lsetvals[0])*x*y

r=0.6

levelset = sqrt(x*x+y*y)-r
n_ref = 7
order = 1
errors = []

for i in range(n_ref):
  V = H1(mesh,order=1)
  lset_approx = GridFunction(V)
  InterpolateToP1(levelset,lset_approx)
  Draw(lset_approx)

  f = CoefficientFunction(1)

  integral = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=f,order=order,domain_type=IF,heapsize=1000000)
  print("Result of Integration Reflevel ",i,": ", integral)
  errors.append(abs(integral - r*pi/2))

  if i < n_ref - 1:
    mesh.Refine()

print("L2-errors:", errors)
l2_eoc = [log(errors[i+1]/errors[i])/log(0.5) for i in range(n_ref-1)]
print("experimental order of convergence (L2):", l2_eoc)
