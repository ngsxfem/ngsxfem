# unfitted convection diffusion equation with Neumann b.c.
from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from numpy import pi

from xfem.lset_spacetime import *

# SetTestoutFile("test.out")
SetNumThreads(1)

square = SplineGeometry()
square.AddRectangle([-0.6,-0.6],[0.6,1.5],bc=1)
ngmesh = square.GenerateMesh(maxh=0.5, quad_dominated=False)
mesh = Mesh (ngmesh)

n_ref = 4
for i in range(n_ref):
    mesh.Refine()
    
h = specialcf.mesh_size

k=2
k_t = k
k_s = k
fes1  = H1(mesh,order=k_s,dirichlet=[])
fes_lset_slice = H1(mesh, order=1, dirichlet=[])
fes_dfm_slice = H1(mesh, order=k_s, dim=mesh.dim)

lset_order_time = k
time_order = 2*k_t
tfe = ScalarTimeFE(k_t) 

st_fes = SpaceTimeFESpace(fes1,tfe,flags = {"dgjumps": True, "no_low_order_space" : True})

tend = 0.5

def SolveProblem(delta_t):
    
    fes1.Update()
    fes_lset_slice.Update()
    fes_dfm_slice.Update()
    st_fes.Update()

    tnew = 0
    
    told = Parameter(0)
    tref = ReferenceTimeVariable()
    t = told + delta_t*tref
    
    lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = lset_order_time,
                                                    threshold=0.05, discontinuous_qn=False,periodic=False)
    
    r0 = 0.5
    r1 = 0.5*pi/r0                                        
    rho =  CoefficientFunction((1/(pi))*sin(2*pi*t))
    d_rho = CoefficientFunction(2*cos(2*pi*t))
    w = CoefficientFunction((0,d_rho)) # convection
    levelset= CoefficientFunction(sqrt(x*x+(y-rho)*(y-rho)) -r0)
    coeff_f = CoefficientFunction(  -(pi/r0)*r1*( sin(r1*sqrt(x*x+(y-rho)*(y-rho)))*sin(r1*sqrt(x*x+(y-rho)*(y-rho))) - cos(r1*sqrt(x*x+(y-rho)*(y-rho)))*cos(r1*sqrt(x*x+(y-rho)*(y-rho))) )  + (pi/r0)*cos(r1*sqrt(x*x+(y-rho)*(y-rho)))*sin(r1*sqrt(x*x+(y-rho)*(y-rho)))*(1/sqrt(x*x+(y-rho)*(y-rho))) )
    
    # for monitoring the error
    l2_infty = []
    l2_dt_infty = []
    l2_grad_infty = []
    u_exact = CoefficientFunction(cos(r1*sqrt(x*x+(y-rho)*(y-rho)))*cos(r1*sqrt(x*x+(y-rho)*(y-rho))) )
    u_dt_exact = CoefficientFunction(2*cos(r1*sqrt(x*x+(y-rho)*(y-rho)))*sin(r1*sqrt(x*x+(y-rho)*(y-rho)))*r1*(y-rho)*d_rho*(1/sqrt(x*x+(y-rho)*(y-rho))) )
    gradu_exact = CoefficientFunction( (-2*r1*sin(r1*sqrt(x*x+(y-rho)*(y-rho)))*cos(r1*sqrt(x*x+(y-rho)*(y-rho)))*(x/sqrt(x*x+(y-rho)*(y-rho))) , 
                                        -2*r1*sin(r1*sqrt(x*x+(y-rho)*(y-rho)))*cos(r1*sqrt(x*x+(y-rho)*(y-rho)))*((y-rho)/sqrt(x*x+(y-rho)*(y-rho))) ) )
   
    u0 = GridFunction(st_fes)
    u0_ic = GridFunction(fes1)
    Pu = GridFunction(fes1)

    u = st_fes.TrialFunction()
    v = st_fes.TestFunction()
    
    lset_top = GridFunction(fes_lset_slice)
    lset_bottom = GridFunction(fes_lset_slice)
    
    dfm_top = GridFunction(fes_dfm_slice)
    dfm_back = GridFunction(fes_dfm_slice)
    dfm_forth = GridFunction(fes_dfm_slice)
    
    t_old = 0  
    
    while tend - t_old > delta_t/2:
        
        dfm_back.vec[:] = dfm_top.vec[:]
        dfm = lset_adap_st.CalcDeformation(levelset,told,t_old,delta_t) 
        dfm_top.vec[:] = dfm.vec[lset_order_time*lset_adap_st.ndof_node : (lset_order_time+1)*lset_adap_st.ndof_node]
        dfm_forth.vec[:] = dfm.vec[0 : lset_adap_st.ndof_node]
        lset_p1 = lset_adap_st.lset_p1    
    
        hasneg_spacetime = lset_adap_st.hasneg_spacetime
        hasneg_spacetime |= lset_adap_st.hasif_spacetime
        active_dofs = GetDofsOfElements(st_fes,hasneg_spacetime)
        
        ba_facets = GetFacetsWithNeighborTypes(mesh,a=hasneg_spacetime,b=lset_adap_st.hasif_spacetime)
        
        lset_neg = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}
        lset_pos = { "levelset" : lset_p1, "domain_type" : POS, "subdivlvl" : 0}
        lset_if = { "levelset" : lset_p1, "domain_type" : IF, "subdivlvl" : 0}
        # Put this into LevelSetMeshAdaptation_Spacetime later on
        lset_bottom.vec[:] = lset_p1.vec[0:fes_lset_slice.ndof]
        lset_neg_bottom = { "levelset" : lset_bottom, "domain_type" : NEG, "subdivlvl" : 0}
        lset_pos_bottom = { "levelset" : lset_bottom, "domain_type" : POS, "subdivlvl" : 0}
        lset_top.vec[:] = lset_p1.vec[lset_order_time*fes_lset_slice.ndof : (lset_order_time+1)*fes_lset_slice.ndof]
        lset_neg_top = { "levelset" : lset_top, "domain_type" : NEG, "subdivlvl" : 0}
        lset_pos_top = { "levelset" : lset_top, "domain_type" : POS, "subdivlvl" : 0}
        
        if t_old == 0:
            mesh.SetDeformation(dfm_forth)
            u0_ic.Set(CoefficientFunction(u_exact)) 
            mesh.UnsetDeformation()
        
        a = RestrictedBilinearForm(st_fes,"a", hasneg_spacetime, ba_facets, check_unused=False)
        
        #a = BilinearForm(st_fes, symmetric = False)

         
        told.Set(t_old)   
        #lset_adap_st.v_p1_st.SetOverrideTime(False)
        lset_adap_st.v_def_st.SetOverrideTime(False)
        
#        a += SymbolicBFI(levelset_domain = lset_neg,
#                         form =
#                           u*v,
#                         time_order=time_order)
        
#        a += SymbolicBFI(levelset_domain = lset_neg,
#                         form =
#                         - u*(dt(v) + dt_vec(dfm)*grad(v) )
#                         + delta_t*grad(u)*grad(v)
#                        - delta_t*u*w*grad(v),
#                         time_order=time_order, definedonelements=hasneg_spacetime)
#        
#        a += SymbolicBFI(levelset_domain = lset_neg_top, form = fix_t(u,1)*fix_t(v,1), definedonelements=hasneg_spacetime)
#        
#        a += SymbolicFacetPatchBFI(form = delta_t*(1+delta_t/h)*0.05*1.0/h*1.0/h*(u-u.Other())*(v-v.Other()),time_order=time_order, skeleton=False, definedonelements=ba_facets) 

        lset_neg1 = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}       
        f = LinearForm(st_fes)
        f += SymbolicLFI(levelset_domain = lset_neg1, form = delta_t*coeff_f*v, time_order=time_order)
#        if t_old == 0:
#            f += SymbolicLFI(levelset_domain = lset_neg_bottom,form = u0_ic*fix_t(v,0) )
#        else:
#            Pu.Set(shifted_eval(u0_ic,dfm_back,dfm_forth))
#            f += SymbolicLFI(levelset_domain = lset_neg_bottom,form = Pu*fix_t(v,0) )
            
        mesh.SetDeformation(dfm)
        a.Assemble()
        f.Assemble()
        mesh.UnsetDeformation()
        
        pre = a.mat.Inverse(active_dofs,"pardiso")
        inv = CGSolver(a.mat,pre,printrates=True,maxsteps=5)
        u0.vec.data =  inv * f.vec
           
        # exploiting the nodal property of the time fe:
        u0_ic.vec[:] = u0.vec[k_t*fes1.ndof : (k_t+1)*fes1.ndof]
        
        # measure l_infinity errors 
        grad_u0 = grad(u0)
        #need to scale, because derivative computed on reference interval
        dt_u0 = (1/delta_t)*dt(u0)
        
        time_quad = [0.1*k for k in range(11)]
        times = [t_old + delta_t * xi for xi in time_quad]
        for ti,xi in zip(times,time_quad):
            told.Set(ti) 
            lset_adap_st.v_def_st.SetTime(xi) 
            st_fes.SetTime(xi) 
            lsetp1 = GridFunction(H1(mesh,order=1))
            InterpolateToP1(levelset,lsetp1)
            lset_ti = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
            mesh.SetDeformation(dfm)
            l2_infty.append(sqrt(Integrate(lset_ti,(u_exact-u0)*(u_exact-u0),mesh,order=2*k_s)))
            l2_dt_infty.append(sqrt(Integrate(lset_ti,(u_dt_exact-dt_u0)*(u_dt_exact-dt_u0),mesh,order=2*k_s)))
            l2_grad_infty.append(sqrt(Integrate(lset_ti,InnerProduct(gradu_exact-grad_u0,gradu_exact-grad_u0),mesh,order=2*k_s)))
            mesh.UnsetDeformation()            
        
        told.Set(t_old) 
        st_fes.SetOverrideTime(False)
        lset_adap_st.v_def_st.SetOverrideTime(False)
        
        t_old = t_old + delta_t
        told.Set(t_old)
        
        mesh.SetDeformation(dfm_top)
        l2error = sqrt(Integrate(lset_neg_top,(u_exact-u0_ic)*(u_exact-u0_ic),mesh,order=2*k_s))
        h1serror = sqrt(Integrate(lset_neg_top,InnerProduct(gradu_exact-grad(u0_ic),gradu_exact-grad(u0_ic)),mesh,order=2*k_s))
        mesh.UnsetDeformation()
        print("t = {0}, l2error = {1}".format(t_old,l2error))
        
        Draw(IfPos(-lset_top,u0_ic,float('nan')),mesh,"u")
        
        
    return l2error,h1serror,max(l2_infty),max(l2_grad_infty),max(l2_dt_infty)

with TaskManager():
    SolveProblem(delta_t = 1/32)

def Study_Conv(where, n_ref = 4,delta_t = 1/32):
    errors = {"L2-u(T)" : [],
              "L2-gradu(T)" : [],
              "L2-u-Infinity" : [],
              "L2-grad(u)-Infinity" : [],
              "L2-Dt(u)-Infinity" : []}
    dt_stepsizes = []
    ref_lvl = 0

    with TaskManager():
        while ref_lvl < n_ref:
            if where == "time":
                dt = tend/2**ref_lvl
                e_l2T,e_h1sT,e_l2infty,e_l2inftygrad,e_l2inftydt = SolveProblem(dt)
                dt_stepsizes.append(dt)
            if where == "space":
                e_l2T,e_h1sT,e_l2infty,e_l2inftygrad,e_l2inftydt = SolveProblem(delta_t)
                if ref_lvl < n_ref -1:
                    mesh.Refine()
            if where == "spacetime":
                dt = tend/2**ref_lvl
                e_l2T,e_h1sT,e_l2infty,e_l2inftygrad,e_l2inftydt = SolveProblem(dt)
                if ref_lvl < n_ref -1:
                    mesh.Refine()      
            errors["L2-u(T)"].append(e_l2T)
            errors["L2-gradu(T)"].append(e_h1sT)
            errors["L2-u-Infinity"].append(e_l2infty)
            errors["L2-grad(u)-Infinity"].append(e_l2inftygrad)
            errors["L2-Dt(u)-Infinity"].append(e_l2inftydt)
            ref_lvl +=1

    print("Studying convergence w.r.t. refinements in: " + where )
    print("")
    
    for key, values in errors.items():
        print("Error: " + key)
        print(values)
        eoc = [ log(values[i-1]/values[i])/log(2) for i in range(1,len(values))]
        print("Eoc: " + key)
        print(eoc)
        print("")
      
  
#Study_Conv(where = "time",n_ref = 5)
#Study_Conv(where = "time",n_ref = 6)
