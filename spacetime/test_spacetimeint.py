from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from numpy import pi

square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
ngmesh = square.GenerateMesh(maxh=0.5, quad_dominated=False)
mesh = Mesh (ngmesh)

fes1 = H1(mesh, order=1, dirichlet=[1,2,3,4])
k_t = 0
tfe = ScalarTimeFE(k_t) 

st_fes = SpaceTimeFESpace(fes1,tfe)
st_fes.SetTime(0.5)

raise Exception("Chirstoph this test is old(?) - lsetp1 is not a space time function")

levelset = (sqrt(x*x+y*y) - 1000.5)
lsetp1 = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lsetp1)
Draw(lsetp1)

lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
lset_pos = { "levelset" : lsetp1, "domain_type" : POS, "subdivlvl" : 0}

t = ReferenceTimeVariable()

u = st_fes.TrialFunction()
v = st_fes.TestFunction()
a = BilinearForm(st_fes)
a += SymbolicBFI(levelset_domain = lset_neg, form = t*1.0*grad(u)*grad(v), time_order=1)
# a += SymbolicBFI(levelset_domain = lset_pos, form = 1.0*u*v)
a.Assemble()

# Draw(t,mesh,"time")

