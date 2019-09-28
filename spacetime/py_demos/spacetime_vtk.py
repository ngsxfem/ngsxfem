# unfitted Heat equation with Neumann b.c.
# solved with a P1-DG-in-time space-time discretization
from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from ngsolve.solvers import *
from xfem import *
from math import pi

from xfem.lset_spacetime import *

#ngsglobals.msg_level = 1

square = SplineGeometry()
square.AddRectangle([-1,-1],[1,1])
ngmesh = square.GenerateMesh(maxh=0.5, quad_dominated=False)
mesh = Mesh (ngmesh)

#### expression for the time variable: 

coef_told = Parameter(0)
coef_delta_t = Parameter(0)
tref = ReferenceTimeVariable()
t = coef_told + coef_delta_t*tref

#### the data: 
# radius of disk (the geometry)
r0 = 0.5


#levelset= x - t
levelset= sqrt((x-0.5*t)**2+y**2) - 0.6

# polynomial order in time
k_t = 1
# polynomial order in space
k_s = 2
# spatial FESpace for solution
fesp1 = H1(mesh, order=1)
# polynomial order in time for level set approximation
lset_order_time = 1
# integration order in time
time_order = 2
# time finite element (nodal!)
tfe = ScalarTimeFE(k_t) 
# space-time finite element space
st_fes = SpaceTimeFESpace(fesp1,tfe, flags = {"dgjumps": True})

#Unfitted heat equation example
tend = 1
delta_t = tend/2
coef_delta_t.Set(delta_t)
tnew = 0
told = 0

lset_p1 = GridFunction(st_fes)

lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = lset_order_time,
                                                threshold=0.1, discontinuous_qn=True)

vtk = SpaceTimeVTKOutput(ma=mesh,coefs=[lset_p1,lset_adap_st.deform],names=["lset_p1","deform"],filename="spacetime_vtk_",
                         subdivision_x=3,subdivision_t=3)
vtk.Do(t_start=coef_told.Get(), t_end=coef_told.Get() + coef_delta_t.Get())

while tend - told > delta_t/2:
    SpaceTimeInterpolateToP1(levelset,tref,lset_p1)
    

    dfm = lset_adap_st.CalcDeformation(levelset,tref) 
    
    vtk.Do(t_start=coef_told.Get(), t_end=coef_told.Get() + coef_delta_t.Get())
    
    told = told + delta_t
    coef_told.Set(told)
