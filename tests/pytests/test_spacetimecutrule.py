import pytest
from ngsolve.meshes import *
from netgen.csg import *
from ngsolve import *
from xfem import *
from math import pi
import netgen.meshing as ngm
from netgen.geom2d import SplineGeometry
from xfem.lset_spacetime import *

tref = ReferenceTimeVariable()

@pytest.mark.parametrize("quad", [True, False])
@pytest.mark.parametrize("integrands", [(tref,0.5,0,1),
                                        (tref**3,0.25,0,3),
                                        ((1-tref)**3,0.25,0,3),
                                        (x,0.5,1,0),
                                        (tref*tref*(x*x+y*y),2/9,2,2)])
def test_spacetime_integrate_no_cut(quad, integrands):
    mesh = MakeStructured2DMesh(quads = quad, nx=1, ny=1)    

    f,ref_value, space_order, time_order = integrands
    
    h1fes = H1(mesh,order=1)
    tfe = ScalarTimeFE(1) 
    fes= SpaceTimeFESpace(h1fes,tfe)
    lset_approx = GridFunction(fes)

    lset_approx.vec[:] = -1

    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : NEG},
                         cf=f, mesh=mesh, order = space_order, time_order=time_order)
    print("Integral: ", integral)
    error = abs(integral - ref_value)
    
    assert error < 5e-15


@pytest.mark.parametrize("quad", [True, False])
@pytest.mark.parametrize("domain", [NEG, POS, IF])
def test_spacetime_integrateX_via_straight_cutted_quad2Dplus1D(domain, quad):
    mesh = MakeStructured2DMesh(quads = quad, nx=1, ny=1)    

    tref = ReferenceTimeVariable()
    
    levelset = lambda t : 1 - 2*x - 2*t
    referencevals = { POS : 1./8, NEG : 1 - 1/8, IF : 1.0/2 }

    h1fes = H1(mesh,order=1)
    lset_approx_h1 = GridFunction(h1fes)
    tfe = ScalarTimeFE(1) 
    fes= SpaceTimeFESpace(h1fes,tfe)
    lset_approx = GridFunction(fes)

    InterpolateToP1(levelset(0),lset_approx_h1)
    lset_approx.vec[0:h1fes.ndof].data = lset_approx_h1.vec
    InterpolateToP1(levelset(1),lset_approx_h1)
    lset_approx.vec[h1fes.ndof:2*h1fes.ndof].data = lset_approx_h1.vec

    print(lset_approx.vec)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = 0, time_order=0)
    print("Integral: ", integral)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15

@pytest.mark.parametrize("pitfal1", [False])
@pytest.mark.parametrize("pitfal2", [False])
@pytest.mark.parametrize("pitfal3", [False])

def test_spacetime_model_spacetime(pitfal1, pitfal2, pitfal3):
    square = SplineGeometry()
    square.AddRectangle([0,0],[1,1],bc=1)
    ngmesh = square.GenerateMesh(maxh=0.05, quad_dominated=False)
    mesh = Mesh (ngmesh)
    
    fes1 = V=H1(mesh, order=1, dirichlet=[1,2,3,4])
    k_t = 1
    tfe = ScalarTimeFE(k_t) 
    
    st_fes = SpaceTimeFESpace(fes1,tfe)
    st_fes_ic = SpaceTimeFESpace(fes1,tfe)
    
    tend = 1.0
    delta_t = 1/32
    
    told = Parameter(0)
    tref = ReferenceTimeVariable()
    t = told + delta_t*tref

    u_exact = lambda t: CoefficientFunction( sin(pi*t)*sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y)  )
    coeff_f = CoefficientFunction( pi*cos(pi*t)*sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y)
                               -2*pi*pi*sin(pi*t)*( cos(pi*x)*cos(pi*x)*sin(pi*y)*sin(pi*y)              
                                                   -sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y)
                                                   +cos(pi*y)*cos(pi*y)*sin(pi*x)*sin(pi*x)
                                                  -sin(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y))) 

    u0 = GridFunction(st_fes)
    u0_ic = GridFunction(fes1)
    u = st_fes.TrialFunction()
    v = st_fes.TestFunction()

    # dummy lset domain to call symboliccutbfi instead of usual symbolicbfi...
    levelset = (sqrt(x*x+y*y) - 1000.5)
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lsetp1)
    lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}

    a = BilinearForm(st_fes,symmetric=False)
    a += SymbolicBFI(levelset_domain = lset_neg, form = delta_t*grad(u)*grad(v), time_order=2)
    a += SymbolicBFI(form = fix_tref(u,0)*fix_tref(v,0) )
    a += SymbolicBFI(levelset_domain = lset_neg, form = dt(u)*v, time_order=2)
    a.Assemble()

    t_old = 0
    u0_ic.Set(u_exact(0))
    if pitfal1:
        u0_ic.Set(u_exact(t))
    
    while tend - t_old > delta_t/2:
        f = LinearForm(st_fes)
        f += SymbolicLFI(levelset_domain = lset_neg, form = delta_t*coeff_f*v, time_order=2)
        f += SymbolicLFI(form = u0_ic*fix_tref(v,0))
        if pitfal2:
            f += SymbolicLFI(form = u0_ic*v)
        f.Assemble()
        
        u0.vec.data = a.mat.Inverse(st_fes.FreeDofs(),"") * f.vec
        
        # exploiting the nodal property of the time fe:
        #u0_ic.vec[:] = u0.vec[0:fes1.ndof]
        u0_ic.vec[:].data = u0.vec[fes1.ndof : 2*fes1.ndof]
    
        t_old = t_old + delta_t
        told.Set(t_old)
        
        l2error = sqrt (Integrate ( (u_exact(t_old) -u0_ic)**2, mesh))
        if pitfal3:
            l2error = sqrt (Integrate ( (u_exact(t) -u0_ic)**2, mesh))
                
        print("t = {0}, l2error = {1}".format(t_old,l2error))
        assert l2error < 5e-3
    assert l2error < 2e-4

def test_spacetime_model_spacetime_caller():
    try:
        test_spacetime_model_spacetime(True, False, False)
    except Exception as e:
        if("TimeVariableCoefficientFunction::Evaluate called with a mere space IR" in str(e)):
            print("Failed properly")
        else:
            print('Unexpected exception raised:', e)
            raise Exception("Wrong kind of failure")
    else:
        raise Exception("No failure at all")
    
    try:
        test_spacetime_model_spacetime(False, True, False)
    except Exception as e:
        if("SpaceTimeFE :: CalcShape called with a mere space IR" in str(e)):
            print("Failed properly")
        else:
            print('Unexpected exception raised:', e)
            raise Exception("Wrong kind of failure")
    else:
        raise Exception("No failure at all")
    
    try:
        test_spacetime_model_spacetime(False, False, True)
    except Exception as e:
        if("TimeVariableCoefficientFunction::Evaluate called with a mere space IR" in str(e)):
            print("Failed properly")
        else:
            print('Unexpected exception raised:', e)
            raise Exception("Wrong kind of failure")
    else:
        raise Exception("No failure at all")

# def test_spacetime_spaceP1_timeCGP1():
#     ngsglobals.msg_level = 1

#     square = SplineGeometry()
#     square.AddRectangle([-1,-1],[1,1])
#     ngmesh = square.GenerateMesh(maxh=0.08, quad_dominated=False)
#     mesh = Mesh (ngmesh)

#     coef_told = Parameter(0)
#     coef_delta_t = Parameter(0)
#     tref = ReferenceTimeVariable()
#     t = coef_told + coef_delta_t*tref

#     r0 = 0.5

#     # position shift of the geometry in time
#     rho =  CoefficientFunction((1/(pi))*sin(2*pi*t))
#     rhoL = lambda t:CoefficientFunction((1/(pi))*sin(2*pi*t))
#     #convection velocity:
#     d_rho = CoefficientFunction(2*cos(2*pi*t))
#     w = CoefficientFunction((0,d_rho)) 

#     # level set
#     r = sqrt(x**2+(y-rho)**2)
#     levelset= r - r0

#     # diffusion coefficient
#     alpha = 1

#     # solution and r.h.s.
#     Q = pi/r0   
#     u_exact = cos(Q*r) * sin(pi*t)
#     u_exactL = lambda t: cos(Q*sqrt(x**2+(y-rhoL(t))**2)) * sin(pi*t)
#     coeff_f = (Q/r * sin(Q*r) + (Q**2) * cos(Q*r)) * sin(pi*t) + pi * cos(Q*r) * cos(pi*t)
#     u_init = u_exact

#     # polynomial order in time
#     k_t = 1
#     # polynomial order in space
#     k_s = 1
#     # spatial FESpace for solution
#     fes1 = H1(mesh, order=k_s)
#     # polynomial order in time for level set approximation
#     lset_order_time = k_t
#     # integration order in time
#     time_order = 2*k_t
#     # time finite element (nodal!)
#     tfe = ScalarTimeFE(k_t) 
#     tfe_i = ScalarTimeFE(k_t, skip_first_node=True) # interior
#     tfe_e = ScalarTimeFE(k_t, only_first_node=True) # exterior (inital values)
#     tfe_t = ScalarTimeFE(k_t-1)                     # test

#     # space-time finite element space
#     st_fes = SpaceTimeFESpace(fes1,tfe, flags = {"dgjumps": True})
#     st_fes_i = SpaceTimeFESpace(fes1,tfe_i, flags = {"dgjumps": True})
#     st_fes_e = SpaceTimeFESpace(fes1,tfe_e, flags = {"dgjumps": True})
#     st_fes_t = SpaceTimeFESpace(fes1,tfe_t, flags = {"dgjumps": True})

#     #Fitted heat equation example
#     tend = 1
#     delta_t = tend/64
#     coef_delta_t.Set(delta_t)
#     tnew = 0
#     told = 0

#     lset_p1 = GridFunction(st_fes)

#     SpaceTimeInterpolateToP1(levelset,tref,lset_p1)

#     lset_top = CreateTimeRestrictedGF(lset_p1,1.0)
#     lset_bottom = CreateTimeRestrictedGF(lset_p1,0.0)

#     gfu_i = GridFunction(st_fes_i)
#     gfu_e = GridFunction(st_fes_e)

#     u_last = CreateTimeRestrictedGF(gfu_e,0)
#     SpaceTimeWeakSet(gfu_e, u_exactL(0.0), fes1)

#     u_i = st_fes_i.TrialFunction()
#     u_e = st_fes_e.TrialFunction()
#     v_t = st_fes_t.TestFunction()

#     h = specialcf.mesh_size

#     lset_neg = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}
#     lset_neg_bottom = { "levelset" : lset_bottom, "domain_type" : NEG, "subdivlvl" : 0}
#     lset_neg_top = { "levelset" : lset_top, "domain_type" : NEG, "subdivlvl" : 0}

#     def SpaceTimeNegBFI(form):
#         return SymbolicBFI(levelset_domain = lset_neg, form = form, time_order=time_order)

#     ci = CutInfo(mesh,time_order=time_order)

        
#     hasneg_integrators_a_i = []
#     hasneg_integrators_a_e = []
#     hasneg_integrators_f = []
#     patch_integrators_a_i = []
#     patch_integrators_a_e = []

#     for hasneg_integrators_a,u in [(hasneg_integrators_a_i,u_i),(hasneg_integrators_a_e,u_e)]:
#         hasneg_integrators_a.append(SpaceTimeNegBFI(form = -dt(v_t)*u))
#         hasneg_integrators_a.append(SpaceTimeNegBFI(form = -delta_t*InnerProduct(w,grad(v_t))*u))
#         hasneg_integrators_a.append(SpaceTimeNegBFI(form = delta_t*alpha*grad(u)*grad(v_t)))

#     for patch_integrators_a,u in [(patch_integrators_a_i,u_i),(patch_integrators_a_e,u_e)]:
#         patch_integrators_a.append(SymbolicFacetPatchBFI(form = delta_t*1.05*h**(-2)*(u-u.Other())*(v_t-v_t.Other()),
#                                                         skeleton=False, time_order=time_order))

#     hasneg_integrators_a_i.append(SymbolicBFI(levelset_domain = lset_neg_top, form = fix_tref(u_i,1)*fix_tref(v_t,1)))

#     hasneg_integrators_a_e.append(SymbolicBFI(levelset_domain = lset_neg_bottom, form = -fix_tref(u_e,0)*fix_tref(v_t,0)))
#     #hasneg_integrators_f.append(SymbolicLFI(levelset_domain = lset_neg_bottom,form = u_last*fix_tref(v,0)))

#     hasneg_integrators_f.append(SymbolicLFI(levelset_domain = lset_neg, form = delta_t*coeff_f*v_t, time_order=time_order)) 

#     a_i = BilinearForm(trialspace = st_fes_i, testspace = st_fes_t, check_unused=False, symmetric=False)
#     for integrator in hasneg_integrators_a_i + patch_integrators_a_i:
#         a_i += integrator

#     a_e = BilinearForm(trialspace = st_fes_e, testspace = st_fes_t, check_unused=False, symmetric=False)
#     for integrator in hasneg_integrators_a_e + patch_integrators_a_e:
#         a_e += integrator

#     f = LinearForm(st_fes_t)

#     for integrator in hasneg_integrators_f:
#         f += integrator

#     while tend - told > delta_t/2:
#         SpaceTimeInterpolateToP1(levelset,tref,lset_p1)
#         RestrictGFInTime(spacetime_gf=lset_p1,reference_time=0.0,space_gf=lset_bottom)
#         RestrictGFInTime(spacetime_gf=lset_p1,reference_time=1.0,space_gf=lset_top)

#         # update markers in (space-time) mesh
#         ci.Update(lset_p1,time_order=time_order)

#         # re-compute the facets for stabilization:
#         ba_facets = GetFacetsWithNeighborTypes(mesh,a=ci.GetElementsOfType(HASNEG),
#                                                     b=ci.GetElementsOfType(IF))
#         # re-evaluate the "active dofs" in the space time slab
#         active_dofs = GetDofsOfElements(st_fes,ci.GetElementsOfType(HASNEG))

#         # re-set definedonelements-markers according to new markings:
#         for integrator in hasneg_integrators_a_i + hasneg_integrators_a_e + hasneg_integrators_f:
#             integrator.SetDefinedOnElements(ci.GetElementsOfType(HASNEG))
#         for integrator in patch_integrators_a:
#             integrator.SetDefinedOnElements(ba_facets)

#         # assemble linear system
#         a_i.Assemble()
#         a_e.Assemble()
#         f.Assemble()

#         # solve linear system
#         inv = a_i.mat.Inverse(active_dofs,inverse="")
#         f.vec.data -= a_e.mat * gfu_e.vec
#         gfu_i.vec.data =  inv * f.vec

#         # evaluate upper trace of solution for
#         #  * for error evaluation 
#         #  * upwind-coupling to next time slab
#         RestrictGFInTime(spacetime_gf=gfu_i,reference_time=1.0,space_gf=u_last)   
        
#         SpaceTimeWeakSet(gfu_e, u_last, fes1)
        
#         # update time variable (float and ParameterCL)
#         told = told + delta_t
#         coef_told.Set(told)
        
#         # compute error at end of time slab
#         l2error = sqrt(Integrate(lset_neg_top,(u_exactL(told) -u_last)**2,mesh))
#         # print time and error
#         print("t = {0:10}, l2error = {1:20}".format(told,l2error),end="\n")
#         assert(l2error < 0.3)
#     assert(l2error < 0.08)

def test_spacetime_spaceP1_timeDGP1():
    ngsglobals.msg_level = 1

    square = SplineGeometry()
    square.AddRectangle([-1,-1],[1,1])
    ngmesh = square.GenerateMesh(maxh=0.08, quad_dominated=False)
    mesh = Mesh (ngmesh)

    coef_told = Parameter(0)
    coef_delta_t = Parameter(0)
    tref = ReferenceTimeVariable()
    t = coef_told + coef_delta_t*tref

    r0 = 0.5

    # position shift of the geometry in time
    rho =  CoefficientFunction((1/(pi))*sin(2*pi*t))
    rhoL = lambda t:CoefficientFunction((1/(pi))*sin(2*pi*t))
    #convection velocity:
    d_rho = CoefficientFunction(2*cos(2*pi*t))
    w = CoefficientFunction((0,d_rho)) 

    # level set
    r = sqrt(x**2+(y-rho)**2)
    levelset= r - r0

    # diffusion coefficient
    alpha = 1

    # solution and r.h.s.
    Q = pi/r0   
    u_exact = cos(Q*r) * sin(pi*t)
    u_exactL = lambda t: cos(Q*sqrt(x**2+(y-rhoL(t))**2)) * sin(pi*t)
    coeff_f = (Q/r * sin(Q*r) + (Q**2) * cos(Q*r)) * sin(pi*t) + pi * cos(Q*r) * cos(pi*t)
    u_init = u_exact

    # polynomial order in time
    k_t = 1
    # polynomial order in space
    k_s = 1
    # spatial FESpace for solution
    fes1 = H1(mesh, order=k_s)
    # polynomial order in time for level set approximation
    lset_order_time = 1
    # integration order in time
    time_order = 2
    # time finite element (nodal!)
    tfe = ScalarTimeFE(k_t) 
    # space-time finite element space
    st_fes = SpaceTimeFESpace(fes1,tfe, flags = {"dgjumps": True})

    #Unfitted heat equation example
    tend = 1
    delta_t = tend/32
    coef_delta_t.Set(delta_t)
    tnew = 0
    told = 0

    lset_p1 = GridFunction(st_fes)

    SpaceTimeInterpolateToP1(levelset,tref,lset_p1)

    lset_top = CreateTimeRestrictedGF(lset_p1,1.0)
    lset_bottom = CreateTimeRestrictedGF(lset_p1,0.0)

    gfu = GridFunction(st_fes)

    u_last = CreateTimeRestrictedGF(gfu,0)
    u_last.Set(u_exactL(0.))

    u,v = st_fes.TnT()

    h = specialcf.mesh_size

    lset_neg = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}
    lset_neg_bottom = { "levelset" : lset_bottom, "domain_type" : NEG, "subdivlvl" : 0}
    lset_neg_top = { "levelset" : lset_top, "domain_type" : NEG, "subdivlvl" : 0}

    def SpaceTimeNegBFI(form):
        return SymbolicBFI(levelset_domain = lset_neg, form = form, time_order=time_order)

    ci = CutInfo(mesh,time_order=time_order)

        
    hasneg_integrators_a = []
    hasneg_integrators_f = []
    patch_integrators_a = []
    hasneg_integrators_a.append(SpaceTimeNegBFI(form = delta_t*alpha*grad(u)*grad(v)))
    hasneg_integrators_a.append(SymbolicBFI(levelset_domain = lset_neg_top, form = fix_tref(u,1)*fix_tref(v,1)))
    hasneg_integrators_a.append(SpaceTimeNegBFI(form = -u*dt(v)))
    hasneg_integrators_a.append(SpaceTimeNegBFI(form = -delta_t*u*InnerProduct(w,grad(v))))
    patch_integrators_a.append(SymbolicFacetPatchBFI(form = delta_t*1.05*h**(-2)*(u-u.Other())*(v-v.Other()),
                                                    skeleton=False, time_order=time_order))
    hasneg_integrators_f.append(SymbolicLFI(levelset_domain = lset_neg, form = delta_t*coeff_f*v, time_order=time_order)) 
    hasneg_integrators_f.append(SymbolicLFI(levelset_domain = lset_neg_bottom,form = u_last*fix_tref(v,0)))


    a = BilinearForm(st_fes,check_unused=False,symmetric=False)
    for integrator in hasneg_integrators_a + patch_integrators_a:
        a += integrator

    f = LinearForm(st_fes)

    for integrator in hasneg_integrators_f:
        f += integrator

    while tend - told > delta_t/2:
        SpaceTimeInterpolateToP1(levelset,tref,lset_p1)
        RestrictGFInTime(spacetime_gf=lset_p1,reference_time=0.0,space_gf=lset_bottom)
        RestrictGFInTime(spacetime_gf=lset_p1,reference_time=1.0,space_gf=lset_top)

        # update markers in (space-time) mesh
        ci.Update(lset_p1,time_order=time_order)

        # re-compute the facets for stabilization:
        ba_facets = GetFacetsWithNeighborTypes(mesh,a=ci.GetElementsOfType(HASNEG),
                                                    b=ci.GetElementsOfType(IF))
        # re-evaluate the "active dofs" in the space time slab
        active_dofs = GetDofsOfElements(st_fes,ci.GetElementsOfType(HASNEG))

        # re-set definedonelements-markers according to new markings:
        for integrator in hasneg_integrators_a + hasneg_integrators_f:
            integrator.SetDefinedOnElements(ci.GetElementsOfType(HASNEG))
        for integrator in patch_integrators_a:
            integrator.SetDefinedOnElements(ba_facets)

        # assemble linear system
        a.Assemble()
        f.Assemble()

        # solve linear system
        inv = a.mat.Inverse(active_dofs,inverse="")
        gfu.vec.data =  inv * f.vec
        

        # evaluate upper trace of solution for
        #  * for error evaluation 
        #  * upwind-coupling to next time slab
        RestrictGFInTime(spacetime_gf=gfu,reference_time=1.0,space_gf=u_last)   

        # update time variable (float and ParameterCL)
        told = told + delta_t
        coef_told.Set(told)
        
        # compute error at end of time slab
        l2error = sqrt(Integrate(lset_neg_top,(u_exactL(told) -u_last)**2,mesh))
        # print time and error
        print("t = {0:10}, l2error = {1:20}".format(told,l2error),end="\n")
        assert(l2error < 0.085)

@pytest.mark.parametrize("i", [3,4,5])
def length_1D_Int_test():
    length = 1
    mesh = Make1DMesh(n=2**(i), mapping= lambda x : 2*length*x-length)
    
    r0 = 0.4
    r = sqrt(x**2)
    
    # level set
    levelset= r - r0
    
    fes1 = H1(mesh, order=1)

    lset_p1 = GridFunction(fes1)
    
    InterpolateToP1(levelset,lset_p1)
    val_vol = Integrate({ "levelset" : lset_p1, "domain_type" : NEG}, CoefficientFunction(1.0), mesh)
    print("Numerical value: ", val_vol)
    print("Error: ", abs(val_vol - 2*r0))
    assert( abs(val_vol - 2*r0) < 1e-10)

def area_of_a_circle_ST_error(n_steps = 8, i=1):
    length = 1
    mesh = Make1DMesh(n=2**(i), mapping= lambda x : 2*length*x-length)

    coef_told = Parameter(0)
    coef_delta_t = Parameter(0)
    tref = ReferenceTimeVariable()
    t = coef_told + coef_delta_t*tref
    
    r0 = 0.9
    r = sqrt(x**2+t**2)
    
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
        SpaceTimeInterpolateToP1(levelset,tref,lset_p1)
    
        val_vol = Integrate({ "levelset" : lset_p1, "domain_type" : NEG}, CoefficientFunction(1.0), mesh, time_order = time_order)
        val_int = Integrate({ "levelset" : lset_p1, "domain_type" : IF}, CoefficientFunction(1.0), mesh, time_order = time_order)
        #print(val_vol, val_int)
        sum_vol += val_vol*delta_t
        sum_int += val_int*delta_t
        
        told = told + delta_t
        coef_told.Set(told)

    print("SUM VOL: ", sum_vol)
    print("VOL: ", pi*r0**2/2)
    vol_err = abs(sum_vol - pi*r0**2/2)
    print("\t\tDIFF: ", vol_err)
    
    print("SUM INT: ", sum_int)
    print("AREA: ", 2*r0)
    int_err = abs(sum_int - 2*r0)
    print("\t\tDIFF: ",int_err)
    return (vol_err, int_err)

def test_spacetime_area_of_a_circle():
    l2errors_vol = []
    for i in range(6):
        (n_steps,i) =  (2**(i+2), i+1)
        (vol_err, int_err) = area_of_a_circle_ST_error(n_steps, i)
        l2errors_vol.append(vol_err)
        assert int_err < 1e-10
    
    print("L2 (VOL): ", l2errors_vol)
    eocs_vol = [log(l2errors_vol[i-1]/l2errors_vol[i])/log(2) for i in range(1,len(l2errors_vol))]
    print("EOCS (VOL): ", eocs_vol)
    avg = sum(eocs_vol)/len(eocs_vol)
    print("Average: ", avg)
    assert avg > 1.9
    
def area_of_a_sphere_ST_error(n_steps = 8, i=1, structured_mesh=False):
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
        SpaceTimeInterpolateToP1(levelset,tref,lset_p1)
    
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

def area_of_a_sphere_ST_error_ho(n_steps = 8, i=3, structured_mesh=False, k = 3):
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
    
    lsetadap = LevelSetMeshAdaptation_Spacetime(mesh, order_space=k, order_time=k, threshold=0.5, discontinuous_qn=True)
    
    time_order = k
    
    tend = r0/2
    delta_t = tend/n_steps
    coef_delta_t.Set(delta_t)
    told = 0
    
    dQ = delta_t * dCut(lsetadap.levelsetp1[INTERVAL], NEG, time_order=time_order, order = k,
                        deformation=lsetadap.deformation[INTERVAL])
    dG = delta_t * dCut(lsetadap.levelsetp1[INTERVAL], IF, time_order=time_order, order = k,
                        deformation=lsetadap.deformation[INTERVAL])
    
    sum_vol = 0
    sum_int = 0
    for i in range(n_steps):
        lsetadap.CalcDeformation(levelset)
    
        val_vol = Integrate( CF(1.)*dQ, mesh)
        val_int = Integrate( CF(1.)*dG, mesh)
        #print(val_vol, val_int)
        sum_vol += val_vol
        sum_int += val_int
        
        told = told + delta_t
        coef_told.Set(told)

    print("SUM VOL: ", sum_vol)
    print("VOL: ", 11/24*pi*r0**3)
    vol_err = abs(sum_vol - 11/24*pi*r0**3)
    print("\t\tDIFF: ", vol_err)
    
    print("SUM INT: ", sum_int)
    print("AREA: ", pi/12*r0**2*(3*sqrt(3) + 2*pi))
    int_err = abs(sum_int - pi/12*r0**2*(3*sqrt(3) + 2*pi))
    print("\t\tDIFF: ",int_err)
    return (vol_err, int_err)

@pytest.mark.parametrize("structured", [True, False])
def test_spacetime_area_of_a_sphere(structured):
    
    l2errors_vol = []
    l2errors_int = []
    for i in range(6):
        (n_steps,i) =  (2**(i+2), i+1)
        (vol_err, int_err) = area_of_a_sphere_ST_error(n_steps, i, structured)
        l2errors_vol.append(vol_err)
        l2errors_int.append(int_err)
    
    print("L2 (VOL): ", l2errors_vol)
    eocs_vol = [log(l2errors_vol[i-1]/l2errors_vol[i])/log(2) for i in range(1,len(l2errors_vol))]
    print("EOCS (VOL): ", eocs_vol)
    avg = sum(eocs_vol)/len(eocs_vol)
    print("Average: ", avg)
    assert avg > 1.9
    
    print("L2 (INT): ", l2errors_int)
    eocs_int = [log(l2errors_int[i-1]/l2errors_int[i])/log(2) for i in range(1,len(l2errors_int))]
    print("EOCS (INT): ", eocs_int)
    avg = sum(eocs_int)/len(eocs_int)
    print("Average: ", avg)
    assert avg > 1.9

@pytest.mark.parametrize("structured", [True, False])
@pytest.mark.parametrize("k", [2,3,4,5,6,7,8])
def test_spacetime_area_of_a_sphere_ho(structured, k):
    
    l2errors_vol = []
    l2errors_int = []
    for i in range(6):
        (n_steps,i) =  (2**(i+1), i+1)
        (vol_err, int_err) = area_of_a_sphere_ST_error_ho(n_steps, i, structured, k)
        l2errors_vol.append(vol_err)
        l2errors_int.append(int_err)
        if vol_err < 1e-9 or int_err < 1e-9:
            break
    
    print("L2 (VOL): ", l2errors_vol)
    eocs_vol = [log(l2errors_vol[i-1]/l2errors_vol[i])/log(2) for i in range(1,len(l2errors_vol))]
    print("EOCS (VOL): ", eocs_vol)
    avg = sum(eocs_vol)/len(eocs_vol)
    print("Average: ", avg)
    assert avg > k+0.8
    
    print("L2 (INT): ", l2errors_int)
    eocs_int = [log(l2errors_int[i-1]/l2errors_int[i])/log(2) for i in range(1,len(l2errors_int))]
    print("EOCS (INT): ", eocs_int)
    avg = sum(eocs_int)/len(eocs_int)
    print("Average: ", avg)
    assert avg > k+0.8

#ngsxfemglobals.do_naive_timeint = True
#ngsxfemglobals.naive_timeint_order = 2
#ngsxfemglobals.naive_timeint_subdivs = 1
#test_spacetime_area_of_a_sphere_ho(False, 2)

def area_of_a_hypersphere_ST_error(n_steps = 64, i=1, structured_mesh= True):
    if structured_mesh:
        length = 1
        mesh = MakeStructured3DMesh(hexes=False,nx=2**(i),ny=2**(i), nz=2**(i),mapping= lambda x,y,z : (2*length*x-length,2*length*y-length, 2*length*z - length))
    else:
        cube = CSGeometry()
        cube.Add (OrthoBrick(Pnt(-1,-1,-1), Pnt(1,1,1)))
        ngmesh = cube.GenerateMesh(maxh=(1/2)**(i-1), quad_dominated=False)
        mesh = Mesh (ngmesh)

    coef_told = Parameter(0)
    coef_delta_t = Parameter(0)
    tref = ReferenceTimeVariable()
    t = coef_told + coef_delta_t*tref
    
    r0 = 0.9
    r = sqrt(x**2+y**2+z**2+t**2)
    
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
        SpaceTimeInterpolateToP1(levelset,tref,lset_p1)
    
        val_vol = Integrate({ "levelset" : lset_p1, "domain_type" : NEG}, CoefficientFunction(1.0), mesh, time_order = time_order)
        val_int = Integrate({ "levelset" : lset_p1, "domain_type" : IF}, CoefficientFunction(1.0), mesh, time_order = time_order)
        #print(val_vol, val_int)
        sum_vol += val_vol*delta_t
        sum_int += val_int*delta_t
        
        told = told + delta_t
        coef_told.Set(told)

    print("SUM VOL: ", sum_vol)
    print("VOL: ", pi**2/4*r0**4)
    vol_err = abs(sum_vol - pi**2/4*r0**4)
    print("\t\tDIFF: ", vol_err)
    
    print("SUM INT: ", sum_int)
    print("AREA: ", 8/3*pi*r0**3)
    int_err = abs(sum_int - 8/3*pi*r0**3)
    print("\t\tDIFF: ", int_err)
    return (vol_err, int_err)

@pytest.mark.parametrize("structured", [True, False])
def test_spacetime_area_of_a_hypersphere(structured):
    l2errors_vol = []
    l2errors_int = []
    for i in range(3):
        (n_steps,i) =  (2**(i+3), i+2)
        (vol_err, int_err) = area_of_a_hypersphere_ST_error(n_steps, i, structured)
        l2errors_vol.append(vol_err)
        l2errors_int.append(int_err)
    
    print("L2 (VOL): ", l2errors_vol)
    eocs_vol = [log(l2errors_vol[i-1]/l2errors_vol[i])/log(2) for i in range(1,len(l2errors_vol))]
    print("EOCS (VOL): ", eocs_vol)
    avg = sum(eocs_vol)/len(eocs_vol)
    print("Average: ", avg)
    assert avg > 1.9
    
    print("L2 (INT): ", l2errors_int)
    eocs_int = [log(l2errors_int[i-1]/l2errors_int[i])/log(2) for i in range(1,len(l2errors_int))]
    print("EOCS (INT): ", eocs_int)
    avg = sum(eocs_int)/len(eocs_int)
    print("Average: ", avg)
    assert avg > 1.9

def test_spacetime_spaceP4_timeDGP4():
    ngsglobals.msg_level = 1

    # -------------------------------- PARAMETERS ---------------------------------
    # DISCRETIZATION PARAMETERS:

    # Parameter for refinement study:
    i = 2
    n_steps = 2**i
    space_refs = i

    # Polynomial order in time
    k_t = 4
    # Polynomial order in space
    k_s = k_t
    # Polynomial order in time for level set approximation
    lset_order_time = k_t
    # Integration order in time
    time_order = 2 * k_t
    # Time stepping parameters
    tstart = 0
    tend = 0.5
    delta_t = (tend - tstart) / n_steps
    maxh = 0.5
    # Ghost-penalty parameter
    gamma = 0.05
    # Map from reference time to physical time
    told = Parameter(tstart)
    t = told + delta_t * tref

    # PROBLEM SETUP:

    # Outer domain:
    rect = SplineGeometry()
    rect.AddRectangle([-0.6, -1], [0.6, 1])

    # Level set geometry
    # Radius of disk (the geometry)
    R = 0.5
    # Position shift of the geometry in time
    rho = (1 / (pi)) * sin(2 * pi * t)
    # Convection velocity:
    w = CoefficientFunction((0, rho.Diff(t)))
    # Level set
    r = sqrt(x**2 + (y - rho)**2)
    levelset = r - R

    # Diffusion coefficient
    alpha = 1
    # Solution
    u_exact = cos(pi * r / R) * sin(pi * t)
    # R.h.s.
    coeff_f = (u_exact.Diff(t)
            - alpha * (u_exact.Diff(x).Diff(x) + u_exact.Diff(y).Diff(y))
            + w[0] * u_exact.Diff(x) + w[1] * u_exact.Diff(y)).Compile()

    # ----------------------------------- MAIN ------------------------------------
    ngmesh = rect.GenerateMesh(maxh=maxh, quad_dominated=False)
    for j in range(space_refs):
        ngmesh.Refine()
    mesh = Mesh(ngmesh)

    # Spatial FESpace for solution
    fes1 = H1(mesh, order=k_s, dgjumps=True)
    # Time finite element (nodal!)
    tfe = ScalarTimeFE(k_t)
    # (Tensor product) space-time finite element space
    st_fes = tfe * fes1

    # Space time version of Levelset Mesh Adapation object. Also offers integrator
    # helper functions that involve the correct mesh deformation
    lsetadap = LevelSetMeshAdaptation_Spacetime(mesh, order_space=k_s,
                                                order_time=lset_order_time,
                                                threshold=0.5,
                                                discontinuous_qn=True)

    gfu = GridFunction(st_fes)
    u_last = CreateTimeRestrictedGF(gfu, 1)

    scene = DrawDC(lsetadap.levelsetp1[TOP], u_last, 0, mesh, "u_last",
                deformation=lsetadap.deformation[TOP])

    u, v = st_fes.TnT()
    h = specialcf.mesh_size

    ba_facets = BitArray(mesh.nfacet)
    ci = CutInfo(mesh, time_order=0)

    dQ = delta_t * dCut(lsetadap.levelsetp1[INTERVAL], NEG, time_order=time_order,
                        deformation=lsetadap.deformation[INTERVAL],
                        definedonelements=ci.GetElementsOfType(HASNEG))
    dOmold = dCut(lsetadap.levelsetp1[BOTTOM], NEG,
                deformation=lsetadap.deformation[BOTTOM],
                definedonelements=ci.GetElementsOfType(HASNEG), tref=0)
    dOmnew = dCut(lsetadap.levelsetp1[TOP], NEG,
                deformation=lsetadap.deformation[TOP],
                definedonelements=ci.GetElementsOfType(HASNEG), tref=1)
    dw = delta_t * dFacetPatch(definedonelements=ba_facets, time_order=time_order,
                            deformation=lsetadap.deformation[INTERVAL])


    def dt(u):
        return 1.0 / delta_t * dtref(u)


    a = RestrictedBilinearForm(st_fes, "a", check_unused=False,
                            element_restriction=ci.GetElementsOfType(HASNEG),
                            facet_restriction=ba_facets)
    a += v * (dt(u) - dt(lsetadap.deform) * grad(u)) * dQ
    a += (alpha * InnerProduct(grad(u), grad(v))) * dQ
    a += (v * InnerProduct(w, grad(u))) * dQ
    a += u * v * dOmold
    a += h**(-2) * (1 + delta_t / h) * gamma * \
        (u - u.Other()) * (v - v.Other()) * dw

    f = LinearForm(st_fes)
    f += coeff_f * v * dQ
    f += u_last * v * dOmold

    # Set initial values
    u_last.Set(fix_tref(u_exact, 0))
    # Project u_last at the beginning of each time step
    lsetadap.ProjectOnUpdate(u_last)

    while tend - told.Get() > delta_t / 2:
        lsetadap.CalcDeformation(levelset)

        # Update markers in (space-time) mesh
        ci.Update(lsetadap.levelsetp1[INTERVAL], time_order=0)

        # re-compute the facets for stabilization:
        ba_facets[:] = GetFacetsWithNeighborTypes(mesh,
                                                a=ci.GetElementsOfType(HASNEG),
                                                b=ci.GetElementsOfType(IF))
        active_dofs = GetDofsOfElements(st_fes, ci.GetElementsOfType(HASNEG))

        a.Assemble(reallocate=True)
        f.Assemble()

        # Solve linear system
        inv = a.mat.Inverse(active_dofs, inverse="")
        gfu.vec.data = inv * f.vec.data

        # Evaluate upper trace of solution for
        #  * for error evaluation
        #  * upwind-coupling to next time slab
        RestrictGFInTime(spacetime_gf=gfu, reference_time=1.0, space_gf=u_last)

        # Compute error at final time
        l2error = sqrt(Integrate((u_exact - u_last)**2 * dOmnew, mesh))

        # Update time variable (ParameterCL)
        told.Set(told.Get() + delta_t)
        print("\rt = {0:12.9f}, L2 error = {1:12.9e}".format(told.Get(), l2error))
    assert(l2error < 1e-2)
    return l2error
