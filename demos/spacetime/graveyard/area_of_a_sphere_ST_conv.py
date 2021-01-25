from ngsolve import *
from ngsolve.meshes import *
from netgen.geom2d import SplineGeometry

from xfem import *
from math import pi

from xfem.lset_spacetime import *

ngsglobals.msg_level = 1

def area_of_a_sphere_ST_error(n_steps = 64, i=3, structured_mesh= True):
    if structured_mesh:
        length = 1
        mesh = MakeStructured2DMesh(quads=False,nx=2**(i),ny=2**(i),mapping= lambda x,y : (2*length*x-length,2*length*y-length))
    else:
        square = SplineGeometry()
        square.AddRectangle([-1,-1],[1,1])
        ngmesh = square.GenerateMesh(maxh=(1/2)**(i-1), quad_dominated=False)
        mesh = Mesh (ngmesh)

    coef_told = Parameter(0)
    coef_delta_t = Parameter(0)
    tref = ReferenceTimeVariable()
    t = coef_told + coef_delta_t*tref
    
    r0 = 0.9
    r = sqrt(x**2+y**2+t**2)
    
    # level set
    levelset= r - r0
    
    time_order = 1
    fes1 = H1(mesh, order=1)
    tfe = ScalarTimeFE(time_order)
    st_fes = SpaceTimeFESpace(fes1,tfe)
    
    tend = 1
    delta_t = tend/n_steps
    coef_delta_t.Set(delta_t)
    told = 0

    lset_p1 = GridFunction(st_fes)
    
    sum_vol = 0
    sum_int = 0
    for i in range(n_steps):
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
    vol_err = abs(sum_vol - 2/3*pi*r0**3)
    print("\t\tDIFF: ", vol_err)
    
    print("SUM INT: ", sum_int)
    print("AREA: ", 0.5*pi**2*r0**2)
    int_err = abs(sum_int - 0.5*pi**2*r0**2)
    print("\t\tDIFF: ",int_err)
    return (vol_err, int_err)

l2errors_vol = []
l2errors_int = []
for i in range(7):
    (n_steps,i) =  (2**(i+2), i+1)
    (vol_err, int_err) = area_of_a_sphere_ST_error(n_steps, i, False)
    l2errors_vol.append(vol_err)
    l2errors_int.append(int_err)

eocs_vol = [log(l2errors_vol[i-1]/l2errors_vol[i])/log(2) for i in range(1,len(l2errors_vol))]
print("EOCS (VOL): ", eocs_vol)
print("Average: ", sum(eocs_vol)/len(eocs_vol))

eocs_int = [log(l2errors_int[i-1]/l2errors_int[i])/log(2) for i in range(1,len(l2errors_int))]
print("EOCS (INT): ", eocs_int)
print("Average: ", sum(eocs_int)/len(eocs_int))
