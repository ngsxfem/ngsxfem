# integration on lset domains

from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *
from netgen.geom2d import SplineGeometry

from netgen.csg import CSGeometry, OrthoBrick, Pnt
cube = OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) ).bc(1)
geom = CSGeometry()
geom.Add (cube)
ngmesh = geom.GenerateMesh(maxh=1, quad_dominated=True)
mesh = Mesh(ngmesh)

def binary_pow(x,a):
    if a == 0:
        return 0*x+1
    elif a == 1:
        return x
    else:
        print("Invalid argument a")

#levelset coefficients
c = [ [ [ 1, -2 ], [-2, 0 ] ] , [[-2, 0], [0,0]] ]
#c = [ [ [ 1.798687, -2.123213 ], [-2.657, 0.8352 ] ] , [[-2.3653, 0.2342356], [-86.4586, -40] ]]

levelset_py = lambda x,y,z: sum( [c[alpha][beta][gamma]*binary_pow(x,alpha)*binary_pow(y,beta)*binary_pow(z,gamma) for alpha in [0,1] for beta in [0,1] for gamma in [0,1]] )
levelset = levelset_py(x,y,z)
referencevals = { POS : 0, NEG : 0 }

V = H1(mesh,order=1)
lset_approx = GridFunction(V)
InterpolateToP1(levelset,lset_approx)

n1 = 1.0/sqrt(InnerProduct(grad(lset_approx),grad(lset_approx))) * grad(lset_approx)
f = CoefficientFunction (1.0)
u = CoefficientFunction ((x,0,0))

domains = [NEG,POS]

for i in range(6):
    if i == 0:
        levelset_l = levelset_py(1,x,y)
        n2 = CoefficientFunction ((1,0,0))
        u_l = CoefficientFunction ((1,0,0))
    elif i == 1:
        levelset_l = levelset_py(0,x,y)
        n2 = CoefficientFunction ((-1,0,0))
        u_l = CoefficientFunction ((0,0,0))
    elif i == 2:
        levelset_l = levelset_py(x,1,y)
        n2 = CoefficientFunction ((0,1,0))
        u_l = u
    elif i == 3:
        levelset_l = levelset_py(x,0,y)
        n2 = CoefficientFunction ((0,-1,0))
        u_l = u
    elif i == 4:
        levelset_l = levelset_py(x,y,1)
        n2 = CoefficientFunction ((0,0,1))
        u_l = u
    elif i == 5:
        levelset_l = levelset_py(x,y,0)
        n2 = CoefficientFunction ((0,0,-1))
        u_l = u

    square = SplineGeometry()
    square.AddRectangle([0,0],[1,1],bc=1)
    mesh_l = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))

    V_l = H1(mesh_l,order=1)
    lset_l_approx = GridFunction(V_l)
    InterpolateToP1(levelset_l,lset_l_approx)
    inte = dict()
    for key in domains:
        inte[key] = NewIntegrateX(lset=lset_l_approx,mesh=mesh_l,cf=u_l*n2,order=32,domain_type=key,heapsize=1000000)
        referencevals[key] += inte[key]
        print("i,key: ",i,key," , inte: ", inte[key])
    print("sum POS+NEG: ", inte[POS]+inte[NEG])

I1 = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=-u*n1,order=32,domain_type=IF,heapsize=1000000)
I2 = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=u*n1,order=32,domain_type=IF,heapsize=1000000)
referencevals[POS] += I1
referencevals[NEG] += I2
print("Interface ints (POS,NEG): ",I1,I2)

print("Referencevals: ", referencevals)

errors = dict()
eoc = dict()

for key in domains:
    errors[key] = []
    eoc[key] = []

for order in range(16):
    for key in domains:
        integral = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=f,order=order,domain_type=key,heapsize=100000000)
        print("Integral on Domain ", key, " : ",integral)
        errors[key].append(abs(integral - referencevals[key]))

Draw(levelset, mesh, "lset")
print("L2 Errors:", errors)
