from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *

from netgen.geom2d import SplineGeometry

# geometry

square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.7, quad_dominated=False))

levelset = sqrt(x*x+y*y) - 0.7
order = 1

lset_approx = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lset_approx)
subdivlvl = 0

# lset_approx = GridFunction(H1(mesh,order=order))
# lset_approx.Set(levelset)
# subdivlvl = 3

# extended FESpace 

VhG = H1(mesh, order=order, dirichlet=[])
gfu = GridFunction(VhG)

# coefficients / parameters: 

n = 1.0/sqrt(InnerProduct(grad(lset_approx),grad(lset_approx))) * grad(lset_approx)
h = specialcf.mesh_size

#tangential projection
def P(u):
   return u - (u*n)*n

# mark only these dofs with support at the interface for the solution
iscut = IsCut(mesh,lset_approx,subdivlvl=subdivlvl)
for i in range(VhG.ndof):
    VhG.FreeDofs()[i] = False
for el in VhG.Elements() :
    if (iscut.vec[el.nr]>0.0):
        for dof in el.dofs:
            VhG.FreeDofs()[dof] = True

# expressions of test and trial functions:

u = VhG.TrialFunction()
v = VhG.TestFunction()

# integration domains (and integration parameter "subdivlvl" and "force_intorder")

lset_if  = { "levelset" : lset_approx, "domain_type" : IF , "subdivlvl" : subdivlvl}

# bilinear forms:

a = BilinearForm(VhG, symmetric = True, flags = { })
a += SymbolicBFI(levelset_domain = lset_if , form = P(grad(u)) * P(grad(v)) + u * v)
a += SymbolicBFI(form = 1.0/h*(iscut * grad(u)*n) * (grad(v)*n))
a.Assemble()

f = LinearForm(VhG)
f += SymbolicLFI(levelset_domain = lset_if, form = sin(x) * v)
f.Assemble();

gfu.vec[:] = 0.0
gfu.vec.data = a.mat.Inverse(VhG.FreeDofs()) * f.vec

nan = CoefficientFunction(float('nan'))
Draw(IfPos(iscut-0.5,gfu,nan),mesh,"u")

# Draw(IsCut(mesh,lset_approx),mesh,"iscut")



