from math import pi
# ngsolve stuff
from ngsolve import *
# visualization stuff
from ngsolve.internal import *
# basic xfem functionality
from xfem import *

from netgen.geom2d import SplineGeometry

# geometry
square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.3, quad_dominated=False))

levelset = sqrt(x*x+y*y) - 0.7
order = 1

lset_approx = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lset_approx)
subdivlvl = 0

# extended FESpace 
VhG = H1(mesh, order=order, dirichlet=[])

# overwrite freedofs of VhG to mark only dofs that are involved in the cut problem
ci = CutInfo(mesh, lset_approx)
ba_IF = ci.GetElementsOfType(IF)
cf_IF = BitArrayCF(ba_IF)
freedofs = VhG.FreeDofs()
freedofs &= GetDofsOfElements(VhG,ba_IF)

gfu = GridFunction(VhG)

# coefficients / parameters: 
n = 1.0/sqrt(InnerProduct(grad(lset_approx),grad(lset_approx))) * grad(lset_approx)
h = specialcf.mesh_size

#tangential projection
def P(u):
   return u - (u*n)*n

# expressions of test and trial functions:
u = VhG.TrialFunction()
v = VhG.TestFunction()

# integration domains (and integration parameter "subdivlvl" and "force_intorder")
lset_if  = { "levelset" : lset_approx, "domain_type" : IF , "subdivlvl" : subdivlvl}

# bilinear forms:

a = BilinearForm(VhG, symmetric = True)
a += SymbolicBFI(levelset_domain = lset_if , form = P(grad(u)) * P(grad(v)) + u * v)
a += SymbolicBFI(form = 1.0/h*( InnerProduct(grad(u),n) * InnerProduct(grad(v),n)), definedonelements = ba_IF)
a.Assemble()

f = LinearForm(VhG)
f += SymbolicLFI(levelset_domain = lset_if, form = sin(x) * v)
f.Assemble();

gfu.vec[:] = 0.0
gfu.vec.data = a.mat.Inverse(freedofs) * f.vec

import sys
if not hasattr(sys, 'argv') or len(sys.argv) == 1 or sys.argv[1] != "testmode":
   nan = CoefficientFunction(float('nan'))
   Draw(IfPos(cf_IF-0.5,gfu,nan),mesh,"u")

   visoptions.mminval = -0.2
   visoptions.mmaxval = 0.2
   visoptions.deformation = 1
   visoptions.autoscale = 0



