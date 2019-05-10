from ngsolve import *
from ngsolve.meshes import *
from netgen.geom2d import SplineGeometry

from xfem import *
from math import pi
from xfem.lsetcurv import *

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
    
    time_order = 5
    space_order = time_order
    
    lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = space_order, order_time = time_order,
                                                threshold=10.5, discontinuous_qn=False)
    
    coef_told = Parameter(0)
    coef_delta_t = Parameter(0)
    tref = ReferenceTimeVariable()
    t = coef_told + coef_delta_t*tref
    
    r0 = 0.9
    r = sqrt(x**2+y**2+t**2)
    
    ref_vals = { IF: pi/12*(3*sqrt(3) + 2*pi)*r0**2, NEG: 11/24*pi*r0**3}
    tend = r0/2
    
    #ref_vals = { IF: 0.5*pi**2*r0**2, NEG: 2/3*pi*r0**3}
    #tend = 1
    
    # level set
    levelset= r - r0
    
    fes1 = H1(mesh, order=space_order)
    tfe = ScalarTimeFE(time_order)
    st_fes = SpaceTimeFESpace(fes1,tfe)
    
    
    delta_t = tend/n_steps
    coef_delta_t.Set(delta_t)
    told = 0

    lset_p1 = lset_adap_st.lset_p1
    
    sum_vol = 0
    sum_int = 0
    for i in range(n_steps):
        #SpaceTimeInterpolateToP1(levelset,tref,0.,delta_t,lset_p1) # call for the master spacetime_weihack -- 0 and tend are ununsed parameter
        
        deformation = lset_adap_st.CalcDeformation(levelset,tref,told, delta_t)
        
        mesh.SetDeformation(deformation)
        #mesh.SetDeformation(test_deformation)
        val_vol = Integrate({ "levelset" : lset_p1, "domain_type" : NEG}, CoefficientFunction(1.0), mesh, order= 2*space_order, time_order = 2*time_order)
        val_int = Integrate({ "levelset" : lset_p1, "domain_type" : IF}, CoefficientFunction(1.0), mesh, order = 2*space_order, time_order = 2*time_order)
        #print(val_vol, val_int)
        mesh.UnsetDeformation()
        
        sum_vol += val_vol*delta_t
        sum_int += val_int*delta_t
        
        told = told + delta_t
        coef_told.Set(told)

    print("SUM VOL: ", sum_vol)
    print("VOL: ", ref_vals[NEG])
    vol_err = abs(sum_vol - ref_vals[NEG])
    print("\t\tDIFF: ", vol_err)
    
    print("SUM INT: ", sum_int)
    print("AREA: ", ref_vals[IF])
    int_err = abs(sum_int - ref_vals[IF])
    print("\t\tDIFF: ",int_err)
    return (vol_err, int_err)

l2errors_vol = []
l2errors_int = []
for i in [0,1,2,3]:
    (n_steps,i) =  (2**(i+2), i+1)
    (vol_err, int_err) = area_of_a_sphere_ST_error(n_steps, i, True)
    l2errors_vol.append(vol_err)
    l2errors_int.append(int_err)

eocs_vol = [log(l2errors_vol[i-1]/l2errors_vol[i])/log(2) for i in range(1,len(l2errors_vol))]
print("EOCS (VOL): ", eocs_vol)
print("Average: ", sum(eocs_vol)/len(eocs_vol))

eocs_int = [log(l2errors_int[i-1]/l2errors_int[i])/log(2) for i in range(1,len(l2errors_int))]
print("EOCS (INT): ", eocs_int)
print("Average: ", sum(eocs_int)/len(eocs_int))
