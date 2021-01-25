from ngsolve import *
from netgen.geom2d import unit_square

from xfem import *
from math import pi
from xfem.lset_spacetime import *

ngmesh = unit_square.GenerateMesh(maxh=0.1, quad_dominated=False)
mesh = Mesh (ngmesh)

coef_told = Parameter(0)
coef_delta_t = Parameter(0)
tref = ReferenceTimeVariable()
t = coef_told + coef_delta_t*tref

k_s = k_t = 2
fes1 = H1(mesh, order=k_s)

tfe = ScalarTimeFE(k_t) 
tfe_i = ScalarTimeFE(k_t, skip_first_node=True) # interior
tfe_e = ScalarTimeFE(k_t, only_first_node=True) # exterior (inital values)

st_fes = SpaceTimeFESpace(fes1,tfe, flags = {"dgjumps": True})
st_fes_i = SpaceTimeFESpace(fes1,tfe_i, flags = {"dgjumps": True})
st_fes_e = SpaceTimeFESpace(fes1,tfe_e, flags = {"dgjumps": True})

gfu = GridFunction(st_fes)
gfu_i = GridFunction(st_fes_i)
gfu_e = GridFunction(st_fes_e)

gfu_to_test = gfu_e

gfu_to_test.Set( t )

for t in [0, 0.5, 1]: #Those are the Gauss-Lobatto points for k=2
    print("t: ", t)
    gfu_slice = CreateTimeRestrictedGF(gfu_to_test,t)
    
    avg = Integrate(gfu_slice,mesh)
    print("avg: ", avg)
