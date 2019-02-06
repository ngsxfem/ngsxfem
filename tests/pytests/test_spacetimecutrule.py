import pytest
from ngsolve import *
from ngsolve.meshes import *
from xfem import *
from math import pi
import netgen.meshing as ngm
from netgen.geom2d import SplineGeometry

tref = ReferenceTimeVariable()

@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("integrands", [(tref,0.5,0,1),
                                        (tref**3,0.25,0,3),
                                        ((1-tref)**3,0.25,0,3),
                                        (x,0.5,1,0),
                                        (tref*tref*(x*x+y*y),2/9,2,2)])
def test_spacetime_integrate_no_cut(quad_dominated, integrands):
    mesh = MakeStructured2DMesh(quads = quad_dominated, nx=1, ny=1)    

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


@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("domain", [NEG, POS, IF])
def test_spacetime_integrateX_via_straight_cutted_quad2Dplus1D(domain, quad_dominated):
    mesh = MakeStructured2DMesh(quads = quad_dominated, nx=1, ny=1)    

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

def test_spacetime_model_spacetime():
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
    a += SymbolicBFI(form = fix_t(u,0)*fix_t(v,0) )
    a += SymbolicBFI(levelset_domain = lset_neg, form = dt(u)*v, time_order=2)
    a.Assemble()

    t_old = 0
    u0_ic.Set(u_exact(0))
    #u0_ic.Set(u_exact(t))
    
    while tend - t_old > delta_t/2:
        f = LinearForm(st_fes)
        f += SymbolicLFI(levelset_domain = lset_neg, form = delta_t*coeff_f*v, time_order=2)
        f += SymbolicLFI(form = u0_ic*fix_t(v,0))
        #f += SymbolicLFI(form = u0_ic*v)
        f.Assemble()
        
        u0.vec.data = a.mat.Inverse(st_fes.FreeDofs(),"umfpack") * f.vec
        
        # exploiting the nodal property of the time fe:
        #u0_ic.vec[:] = u0.vec[0:fes1.ndof]
        u0_ic.vec[:].data = u0.vec[fes1.ndof : 2*fes1.ndof]
    
        t_old = t_old + delta_t
        told.Set(t_old)
        
        l2error = sqrt (Integrate ( (u_exact(t_old) -u0_ic)**2, mesh))
        #l2error = sqrt (Integrate ( (u_exact(t) -u0_ic)**2, mesh))
                
        print("t = {0}, l2error = {1}".format(t_old,l2error))
        assert l2error < 5e-3
    assert l2error < 2e-4
