# integration on lset domains

from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *
from sympy import *

from netgen.geom2d import SplineGeometry

square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
#mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=False))
mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))

#levelset values on the 4 vertices of the quad
lsetvals = [-0.5832,-1.2,1.234521,0.89427] #case b)
#lsetvals = [-0.18687,0.765764,0.324987,0.48983]
#lsetvals = [0.765764,0.324987,0.48983, -0.18687]
#lsetvals = [1,2/3,-1,-2/3] #case c)
#lsetvals = [1,-1,1/3,-1] #case d)
d = lsetvals[0]
c = lsetvals[2]-lsetvals[1]-lsetvals[3]+lsetvals[0]
a = lsetvals[1] - lsetvals[0]
b = lsetvals[3] - lsetvals[0]
#a = -2
#b = -2
#c = 0.4
#d = 1
levelset = d+a*x + b*y + c*x*y
f = x+0*y

xs = Symbol('xs')
ys = Symbol('ys')
levelset_py = d+a*xs +b*ys + c*xs*ys
f_py = 1*xs+0*ys

referencevals = {} #{ POS : 1./8, NEG : 7./8, IF : 1./sqrt(2)}

y_ast = -(a*xs+d)/(b+c*xs)
I1 = integrate(integrate(f_py, (ys, 0, y_ast)), (xs,0,1))
I2 = integrate(integrate(f_py, (ys, y_ast, 1)), (xs,0,1))
print("I1:", I1)
print("I2:", I2)

print("Sum of I1 and I2:", I1+I2)

if(lsetvals[0] > 0):
    referencevals[POS] = I1
    referencevals[NEG] = I2
else:
    referencevals[POS] = I2
    referencevals[NEG] = I1

V = H1(mesh,order=1)
lset_approx = GridFunction(V)
InterpolateToP1(levelset,lset_approx)

domains = [NEG,POS]

errors = dict()
eoc = dict()

for key in domains:
    errors[key] = []
    eoc[key] = []
inte = dict()

for order in range(16):
    for key in domains:
        integral = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=f,order=order,domain_type=key,heapsize=1000000)
        inte[key] = integral
        print("Integral on Domain ", key, " : ",integral)
        errors[key].append(abs(integral - referencevals[key]))
    print("Sum of Part NEG, POS: ", inte[NEG]+inte[POS])

Draw(levelset, mesh, "lset")
print("L2 Errors:", errors)
