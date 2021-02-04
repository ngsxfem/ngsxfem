"""
unfitted Heat equation with Neumann b.c. solved with a P1-DG-in-time
space-time discretisation
"""
# ------------------------------ LOAD LIBRARIES -------------------------------
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
# -------------------------------- PARAMETERS ---------------------------------


# ----------------------------------- MAIN ------------------------------------

square = SplineGeometry()
A = 1.25
square.AddRectangle([-A,-A],[A,A])
ngmesh = square.GenerateMesh(maxh=0.3, quad_dominated=False)
mesh = Mesh (ngmesh)

#### expression for the time variable: 

coef_told = Parameter(0)
coef_delta_t = Parameter(0)
tref = ReferenceTimeVariable()
t = coef_told + coef_delta_t*tref

#### the data: 
# radius of disk (the geometry)
r0 = 0.5
x0 = lambda t: r0 * cos(pi*t)
y0 = lambda t: r0 * sin(pi*t)

#levelset= x - t
levelset= sqrt((x-x0(t))**2+(y-y0(t))**2) - 0.4

# polynomial order in time
k_t = 2
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
delta_t = tend/4
coef_delta_t.Set(delta_t)
tnew = 0
told = 0

lset_p1 = GridFunction(st_fes)

lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = lset_order_time,
                                                threshold=0.1, discontinuous_qn=True)

ci = CutInfo(mesh,time_order=0)

vtk = SpaceTimeVTKOutput(ma=mesh,coefs=[levelset,lset_p1,
                                        CoefficientFunction((lset_adap_st.deform[0],lset_adap_st.deform[1],0)),
                                        BitArrayCF(ci.GetElementsOfType(IF))],
                         names=["levelset","lset_p1","deform","cutelements"],
                         filename="spacetime_vtk_",
                         subdivision_x=3,subdivision_t=3)
# vtk.Do(t_start=coef_told.Get(), t_end=coef_told.Get() + coef_delta_t.Get())

while tend - told > delta_t/2:
    SpaceTimeInterpolateToP1(levelset,tref,lset_p1)
    dfm = lset_adap_st.CalcDeformation(levelset,tref)
    ci.Update(lset_p1,time_order=0)
    
    vtk.Do(t_start=coef_told.Get(), t_end=coef_told.Get() + coef_delta_t.Get())
    told = told + delta_t
    coef_told.Set(told)
