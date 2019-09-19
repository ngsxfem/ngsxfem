from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from numpy import pi
from time import sleep

from xfem.lset_spacetime import *


square = SplineGeometry()
square.AddRectangle([-0.6,-0.6],[0.6,1],bc=1)
ngmesh = square.GenerateMesh(maxh=0.5, quad_dominated=False)
mesh = Mesh (ngmesh)

for i in range(0):
    mesh.Refine()

k_s = 2
lset_order_time = 3
fes_lset_slice = H1(mesh, order=k_s, dirichlet=[])

intermediate_times = [0.1*i for i in range(11)]
#intermediate_times = [0.5]
tend = 0.5

def SolveProblem(delta_t):
    
    told = Parameter(0)
    tref = ReferenceTimeVariable()
    t = told + delta_t*tref
    t_old = 0
    
    lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = lset_order_time,
                                             threshold=0.1, discontinuous_qn=True)
# Non smooth version
#    r0 = 0.5
#    rho =  CoefficientFunction((1/(pi))*sin(2*pi*t))
#    levelset= CoefficientFunction(sqrt(x*x+(y-rho)*(y-rho)) -r0)

# smooth version
#    r0 = 0.25
#    rho =  CoefficientFunction((1/(pi))*sin(2*pi*t))
#    levelset= CoefficientFunction(x*x+(y-rho)*(y-rho) -r0)

    r0 = 0.25
    rho =  CoefficientFunction((1/(pi))*sin(2*pi*t))
    levelset= CoefficientFunction(x*x*x+(y-rho)*(y-rho) -r0)
    
    
    l2_error = 0
    
    lset_slice = GridFunction(fes_lset_slice)
    
    while tend - t_old > delta_t/2:
               
        dfm = lset_adap_st.CalcDeformation(levelset,told,t_old,delta_t) 
                    
        #
        for atime in intermediate_times:
        #for i,atime in enumerate(lset_adap_st.v_ho_st.TimeFE_nodes().NumPy()): 
            lset_ho = RestrictToTime(lset_adap_st.lset_ho,atime)
            #lset_slice.vec[:] = lset_adap_st.lset_ho.vec[i*lset_adap_st.ndof_node : (i+1)*lset_adap_st.ndof_node]
            

            #Draw(lset_ho,mesh,"restrictLset")
            told.Set(t_old+delta_t*atime)

           
            l2_error = max(l2_error,sqrt(Integrate((lset_ho-levelset)*(lset_ho-levelset),mesh,order=2*k_s)))
            
            #lset_slice.Set(levelset)
            #l2_error = max(l2_error,sqrt(Integrate((lset_slice-levelset)*(lset_slice-levelset),mesh,order=2*k_s)))
            
            #lset_slice.Set(x*x*1-t)
            #l2_error = max(l2_error,sqrt(Integrate((lset_slice-(x*x*1-t))*(lset_slice-(x*x*1-t)),mesh,order=2*k_s)))
            print(l2_error)       
        
        t_old = t_old + delta_t
        told.Set(t_old)
        
    return l2_error
        
#print("Maximum Error = {0}".format(SolveProblem(delta_t = 1/16)))

    
def Study_Conv(where, n_ref = 5,delta_t = 1/32):
    l2_errors = []
    dof_numbers = []
    dt_stepsizes = []
    ref_lvl = 0

    with TaskManager():
        while ref_lvl < n_ref:
            if where == "time":
                dt = tend/2**ref_lvl
                l2_errors.append(SolveProblem(dt))
                dt_stepsizes.append(dt)
            if where == "space":
                l2_errors.append(SolveProblem(delta_t))
                #dof_numbers.append(fes1.ndof)
                if ref_lvl < n_ref -1:
                    mesh.Refine()
            if where == "spacetime":
                dt = tend/2**ref_lvl
                l2_errors.append(SolveProblem(dt))
                if ref_lvl < n_ref -1:
                    mesh.Refine()                
            ref_lvl +=1

            print("Studying convergence w.r.t. refinements in: " + where )
            print("L2-Errors:" )
            print(l2_errors)
            print("Eoc:")
            eoc = [ log(l2_errors[i-1]/l2_errors[i])/log(2) for i in range(1,len(l2_errors))]
            print(eoc)

Study_Conv(where = "space",n_ref = 5,delta_t = 1/128)
#Study_Conv(where = "spacetime",n_ref = 5)
#Study_Conv(where = "time",n_ref = 5)

