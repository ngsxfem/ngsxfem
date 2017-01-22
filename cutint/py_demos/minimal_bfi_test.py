
from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *

def get_levelset(lsetvals):
    return lsetvals[0] +(lsetvals[1] - lsetvals[0])*x + (lsetvals[3] - lsetvals[0])*y + (lsetvals[2]-lsetvals[1]-lsetvals[3]+lsetvals[0])*x*y

from netgen.geom2d import SplineGeometry
square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))

lsetvals = [1,2,3,-4]
levelset = get_levelset(lsetvals)

lset_approx = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lset_approx)

bulkfes = H1(mesh, order=1, dirichlet=[1,2,3,4])
VhG = FESpace([bulkfes,bulkfes])

Draw(levelset, mesh, "lset")
Draw(lset_approx)

lset_neg = { "levelset" : lset_approx, "domain_type" : NEG, "subdivlvl" : 0}
lset_pos = { "levelset" : lset_approx, "domain_type" : POS, "subdivlvl" : 0}
lset_if  = { "levelset" : lset_approx, "domain_type" : IF , "subdivlvl" : 0}

# bilinear forms:

a = BilinearForm(VhG, symmetric = True, flags = { })
a += SymbolicBFI(levelset_domain = lset_neg, form = x)
a += SymbolicBFI(levelset_domain = lset_pos, form = x)
a += SymbolicBFI(levelset_domain = lset_if , form = x)

a.Assemble()
