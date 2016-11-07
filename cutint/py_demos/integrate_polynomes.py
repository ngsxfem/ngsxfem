from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *

def pow(a,b):
    return exp(log(a)*b)

def Make2DProblem(maxh):
    from netgen.geom2d import SplineGeometry
    square = SplineGeometry()
    square.AddRectangle([0,0],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=maxh, quad_dominated=False))
    return mesh;

def Make3DProblem(maxh):
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    cube = OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) ).bc(1)
    geom = CSGeometry()
    geom.Add (cube)
    mesh = Mesh (geom.GenerateMesh(maxh=maxh, quad_dominated=False))
    return mesh

def MakeProblem(D, maxh):
    if D == 2:
        return Make2DProblem(maxh)
    elif D == 3:
        return Make3DProblem(maxh)

x_bnd = 0.5363452

levelset = x- x_bnd

referencevals = { POS : lambda a,b: 1./((1+a)*(1+b))*(1 - pow(x_bnd, a+1)), NEG : lambda a,b: 1./((1+a)*(1+b))*pow(x_bnd, a+1), IF : lambda a,b: pow(x_bnd, a) }

mesh = MakeProblem(D=2, maxh=0.22)

V = H1(mesh,order=1)
lset_approx = GridFunction(V)
InterpolateToP1(levelset,lset_approx)

domains = [NEG,POS,IF]

errors = dict()
eoc = dict()

for key in domains:
    errors[key] = []
    eoc[key] = []

exponents = [(0,0)] #[(0,1),(1,0),(2,1),(3,3),(6,2)]

for (a,b) in exponents:
    #f = CoefficientFunction(1)
    #f = pow(x,a)*pow(y,b)
    f = 0*x+0*y+1
    
    for key in domains:
        integral = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=f,order=0,domain_type=key,heapsize=1000000)
        errors[key].append(abs(integral - referencevals[key](a,b)))

print("errors:  \n{}\n".format(  errors))
print("   eoc:  \n{}\n".format(     eoc))
