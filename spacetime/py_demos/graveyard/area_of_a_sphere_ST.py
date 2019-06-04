from ngsolve import *
from netgen.geom2d import SplineGeometry

from xfem import *
from math import pi

from xfem.lset_spacetime import *

ngsglobals.msg_level = 1

square = SplineGeometry()
square.AddRectangle([-1,-1],[1,1])
ngmesh = square.GenerateMesh(maxh=0.02, quad_dominated=False)
mesh = Mesh (ngmesh)

coef_told = Parameter(0)
coef_delta_t = Parameter(0)
tref = ReferenceTimeVariable()
t = coef_told + coef_delta_t*tref

r0 = 0.5
r = sqrt(x**2+y**2+t**2)

# level set
levelset= r - r0

time_order = 1
fes1 = H1(mesh, order=1)
tfe = ScalarTimeFE(time_order)
st_fes = SpaceTimeFESpace(fes1,tfe)

n_steps = 100

tend = 1
delta_t = tend/n_steps
coef_delta_t.Set(delta_t)
told = 0

lset_p1 = GridFunction(st_fes)

sum_vol = 0
sum_int = 0
for i in range(n_steps):
    #SpaceTimeInterpolateToP1(levelset,coef_told,told,delta_t,lset_p1) #call for the master Branch
    SpaceTimeInterpolateToP1(levelset,tref,0.,delta_t,lset_p1) # call for the master spacetime_weihack -- 0 and tend are ununsed parameter
    
    val_vol = Integrate({ "levelset" : lset_p1, "domain_type" : NEG}, CoefficientFunction(1.0), mesh, time_order = time_order)
    val_int = Integrate({ "levelset" : lset_p1, "domain_type" : IF}, CoefficientFunction(1.0), mesh, time_order = time_order)
    #print(val_vol, val_int)
    sum_vol += val_vol*delta_t
    sum_int += val_int*delta_t
    
    told = told + delta_t
    coef_told.Set(told)

print("SUM VOL: ", sum_vol)
print("VOL: ", 2/3*pi*r0**3)
print("\t\tDIFF: ", abs(sum_vol - 2/3*pi*r0**3))

print("SUM INT: ", sum_int)
print("AREA: ", 0.5*pi**2*r0**2)
print("\t\tDIFF: ", abs(sum_int - 0.5*pi**2*r0**2))
