import pytest
from ngsolve.meshes import *
from ngsolve import *
from xfem import *
from math import pi
from xfem.lsetcurv import *

@pytest.mark.parametrize("domain", [NEG, POS, IF])
@pytest.mark.parametrize("alpha", [2,4,8])

def test_polynomial_ET_Segm(domain, alpha):
    order = alpha
    mesh = Make1DMesh(1)
    
    x_ast = 0.78522
    levelset = x_ast-x
    referencevals = {POS: x_ast**(alpha+1)/(alpha+1), NEG:(1-x_ast**(alpha+1))/(alpha+1), IF: x_ast**(alpha)}
    V = H1(mesh, order=1)
    lset_approx = GridFunction(V)
    #InterpolateToP1(levelset,lset_approx)
    lset_approx.Set(levelset)
    
    f = x**alpha
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    print("Result of Integration Key ",domain," : ", integral)
    error = abs(integral - referencevals[domain])
    print("Error: ", error)
    
    assert error < 5e-15*(order+1)*(order+1)

@pytest.mark.parametrize("quad", [True, False])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_via_circle_geom(quad, order, domain):
    r=0.6

    levelset = sqrt(x*x+y*y)-r
    referencevals = { POS : 1-pi*r*r/4, NEG : pi*r*r/4, IF : r*pi/2}

    n_ref = 8
    errors = []

    for i in range(n_ref):
        mesh = MakeStructured2DMesh(quads = quad, nx=2**i, ny=2**i)    

        V = H1(mesh,order=1)
        lset_approx = GridFunction(V)
        InterpolateToP1(levelset,lset_approx)
    
        f = CoefficientFunction(1)
    
        integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                             cf=f, mesh=mesh, order = order)
        print("Result of Integration Reflevel ",i,", Key ",domain," : ", integral)
        errors.append(abs(integral - referencevals[domain]))

        
    eoc = [log(errors[i+1]/errors[i])/log(0.5) for i in range(n_ref-1)]

    print("L2-errors:", errors)
    print("experimental order of convergence (L2):", eoc)

    mean_eoc_array = eoc[1:]
    mean_eoc = sum(mean_eoc_array)/len(mean_eoc_array)
    assert mean_eoc > 1.75

@pytest.mark.parametrize("order", [2,4])
@pytest.mark.parametrize("domain", [POS, NEG])
def test_new_integrateX_via_sphere_geom_quad(order, domain):
    r=0.7234436998

    levelset = sqrt(x*x+y*y+z*z)-r
    referencevals = { POS : 1-pi*r*r*r/6, NEG : pi*r*r*r/6, IF : r*r*pi/2}

    n_ref = 6
    errors = []

    for i in range(n_ref):
        mesh = MakeStructured3DMesh(hexes = True, nx=2**i, ny=2**i, nz=2**i)
        print("i: " +str(i))
        print("Argument Meshing: ",str(2**i))
        
        V = H1(mesh,order=1)
        lset_approx = GridFunction(V)
        InterpolateToP1(levelset,lset_approx)
    
        f = CoefficientFunction(1)
    
        integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
        print("Result of Integration Reflevel ",i,", Key ",domain," : ", integral)
        errors.append(abs(integral - referencevals[domain]))
        
    eoc = [log(errors[i+1]/errors[i])/log(0.5) for i in range(n_ref-1)]

    print("L2-errors:", errors)
    print("experimental order of convergence (L2):", eoc)

    mean_eoc_array = eoc[1:]
    mean_eoc = sum(mean_eoc_array)/len(mean_eoc_array)
    assert mean_eoc > 1.75
    
@pytest.mark.parametrize("quad", [True, False])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])
@pytest.mark.parametrize("N", [1,10,30])

def test_new_integrateX_via_straight_cutted_quad2D(order, domain, quad, N):
    mesh = MakeStructured2DMesh(quads = quad, nx=N, ny=N)    
    
    levelset = 1 - 2*x - 2*y
    referencevals = {NEG: 7/8, POS: 1/8, IF: 1/sqrt(2)}
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)

@pytest.mark.parametrize("quad", [False, True])
@pytest.mark.parametrize("order", [2])
@pytest.mark.parametrize("domain", [IF, NEG, POS])
@pytest.mark.parametrize("dim", [x,y])
@pytest.mark.parametrize("eps", [1e-1, 1e-2, 5e-3, 1e-3, 0])

def test_new_integrateX_via_orth_cutted_quad2D_epsiloned(order, domain, quad, dim, eps):
    mesh = MakeStructured2DMesh(quads = quad, nx=1, ny=1)    
    
    if dim == x:
        levelset = 1 - 2*x + eps*(y-0.5)
    elif dim == y:
        levelset = 1 - 2*y + eps*(x-0.5)
        
    referencevals = {NEG: 1./2, POS: 1./2, IF: sqrt(1.+eps*eps/4) }
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    print(error)
    
    assert error < 5e-15*(order+1)*(order+1)


@pytest.mark.parametrize("quad", [True, False])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS])
@pytest.mark.parametrize("alpha", [0,1,2])
@pytest.mark.parametrize("dim", [x,y])

#integrate f(x) = dim^alpha on the geometry implied by phi(x,y,z) = 1 - 2*x - 2*y
# for analytic solution see
# http://www.wolframalpha.com/input/?i=integrate+from+0+to+1%2F2+from+0+to+(1%2F2-x)+x%5Ealpha+dy+dx
def test_new_integrateX_via_straight_cutted_quad2D_polynomial(order, domain, quad, alpha, dim):
    mesh = MakeStructured2DMesh(quads = quad, nx=1, ny=1)    
    # square = SplineGeometry()
    # square.AddRectangle([0,0],[1,1],bc=1)
    # mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=quad))
    
    levelset = 1 - 2*x - 2*y
    val_pos = 2**(-alpha-2)/(alpha*alpha + 3*alpha+2)
    referencevals = {POS: val_pos, NEG: 1./(alpha+1) - val_pos}
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = dim**alpha
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)

@pytest.mark.parametrize("quad", [True, False])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_via_straight_cutted_quad3D(order, domain, quad):
    mesh = MakeStructured3DMesh(hexes = quad, nx=1, ny=1, nz=1)    
    
    levelset = 1 - 2*x - 2*y - 2*z
    referencevals = { POS : 1./48, NEG : 47./48, IF : sqrt(3)/8 }
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    print("Integral: ", integral)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)

@pytest.mark.parametrize("quad", [True, False])
@pytest.mark.parametrize("order", [4])
@pytest.mark.parametrize("domain", [NEG, POS])
@pytest.mark.parametrize("alpha", [0,1,2])
@pytest.mark.parametrize("dim", [x,y,z])

#integrate f(x) = dim^alpha on the geometry implied by phi(x,y,z) = 1 - 2*x - 2*y - 2*z
# for analytic solution see
# http://www.wolframalpha.com/input/?i=integrate+from+0+to+1%2F2+from+0+to+(1%2F2-x)+from+0+to+(1%2F2-x-y)+x%5Ealpha+dz+dy+dx
def test_new_integrateX_via_straight_cutted_quad3D_polynomial(order, domain, quad, alpha, dim):
    ngsglobals.msg_level = 0
    mesh = MakeStructured3DMesh(hexes = quad, nx=1, ny=1, nz=1)    
        
    levelset = 1 - 2*x- 2*y - 2*z
    val_pos = 2**(-alpha-3)/(alpha**3+6*alpha*alpha + 11*alpha+6)
    referencevals = {POS: val_pos, NEG: 1./(alpha+1) - val_pos}
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = dim**alpha
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)

@pytest.mark.parametrize("quad", [True, False])
@pytest.mark.parametrize("order", [2])
@pytest.mark.parametrize("domain", [NEG, POS, IF])
@pytest.mark.parametrize("dim", [x,y])

def test_new_integrateX_via_orth_cutted_quad2D(order, domain, quad, dim):
    mesh = MakeStructured2DMesh(quads = quad, nx=1, ny=1)    
    
    levelset = 1 - 3*dim
    referencevals = {NEG: 2./3, POS: 1./3, IF: 1. }
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)
    

@pytest.mark.parametrize("quad", [False, True])
@pytest.mark.parametrize("order", [2])
@pytest.mark.parametrize("domain", [NEG, POS])
@pytest.mark.parametrize("dim", [x,y,z])

def test_new_integrateX_via_orth_cutted_quad3D(order, domain, quad, dim):
    mesh = MakeStructured3DMesh(hexes = quad, nx=1, ny=1, nz=1)    
    
    levelset = 1 - 2*dim
    referencevals = { POS : 1./2, NEG : 1./2, IF : 1. }
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)

def test_eb_cut_integrator_2d():
    from netgen.geom2d import SplineGeometry

    ngsglobals.msg_level = 1

    levelset = (x -1) **4 + (y -1)**4 + 12*y + 1.4*(x -2)**3 - 15
    len_box = 4

    exact = sin(y)
    exact_grad = CoefficientFunction((0, cos(y)))

    order = 5

    l2errors = []
    h1errors = []

    extend_Fh = False
    condense = True

    for i in [2,3,4]:
        square = SplineGeometry()
        square.AddRectangle([-len_box,-len_box],[len_box,len_box],bc=1)
        mesh = Mesh (square.GenerateMesh(maxh=len_box/1.5*0.5**(i), quad_dominated=False))
        print("Max_h = ", len_box/1.5*0.5**(i))

        lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=10.25, discontinuous_qn=True)
        deformation = lsetmeshadap.CalcDeformation(levelset)
        lsetp1 = lsetmeshadap.lset_p1
            
        # TraceFESpace 
        fes1 = L2(mesh, order=order)
        fes2 = FacetFESpace(mesh, order=order)
        fes3 = FacetFESpace(mesh, order=order - 1)
        
        Vhg = FESpace([fes1, fes2, fes3], dgjumps=not condense)
        
        # overwrite freedofs of VhG to mark only dofs that are involved in the cut problem
        ci = CutInfo(mesh, lsetp1)
        reg_Th = ci.GetElementsOfType(IF)   
        reg_Fh = GetFacetsWithNeighborTypes(mesh,a=reg_Th,b=reg_Th,use_and=not extend_Fh)

        gf_reg_Fh = GridFunction(FacetFESpace(mesh,order=0))
        for i in range(len(reg_Fh)):
            gf_reg_Fh.vec[i] = 1 if reg_Fh[i] else 0
        
        if not extend_Fh:
                if not condense:
                    freedofs = CompoundBitArray( [GetDofsOfElements(fes1,reg_Th),  GetDofsOfFacets(fes2,reg_Fh) ,  GetDofsOfFacets(fes3,reg_Fh) ] )
                else:
                    freedofs = CompoundBitArray( [fes1.FreeDofs(True) & GetDofsOfElements(fes1,reg_Th),  GetDofsOfFacets(fes2,reg_Fh) ,  GetDofsOfFacets(fes3,reg_Fh) ] )
        else:
                if not condense:
                    freedofs = GetDofsOfElements(Vhg,reg_Th)
                else:
                    freedofs = CompoundBitArray( [fes1.FreeDofs(True) & GetDofsOfElements(fes1,reg_Th),  GetDofsOfElements(fes2,reg_Th) ,  GetDofsOfElements(fes3,reg_Th) ] )
        
        normal_helper_gf = GridFunction(HDiv(mesh, order=0))
        for i in range(len(normal_helper_gf.vec)):
            normal_helper_gf.vec[i] = 1
        
        gfu = GridFunction(Vhg)
            
        #tangential projection to given normal
        def P(u,n_phi):
            return u - (u*n_phi)*n_phi
        
        #normalization (pointwise) of a vector
        def Normalized(u):
            return 1.0 / Norm(u) * u
        
        n_phi1 = Normalized(grad(lsetp1))
        n_phi2 = Normalized(grad(lsetp1).Other())
        
        h = specialcf.mesh_size
        n_F = specialcf.normal(2)
        
        conormal1 = Normalized(P(n_F,n_phi1))
        conormal2 = Normalized(P(-n_F,n_phi2))
        
        normal_helper_proj = InnerProduct(normal_helper_gf, n_F)*n_F
        normal_helper = Normalized(normal_helper_proj)
            
        def jump(u, uhat):
            return u - uhat
        
        beta_E = 4 * (order+1)**2
        beta_F = 100.
        beta_F2 = 1.
        lam_nd = 0.0 if order == 1 else 0.1*1./h + 0.1*h
        
        u, uhat, sigmahat = Vhg.TrialFunction()
        v, vhat, tauhat = Vhg.TestFunction()
        lset_if  = { "levelset" : lsetp1, "domain_type" : IF , "subdivlvl" : 0}
        
        a = RestrictedBilinearForm(Vhg,"a",reg_Th,reg_Fh,check_unused=False, flags={"eliminate_internal": condense})
        a += SymbolicBFI(levelset_domain = lset_if, form = P(grad(u),n_phi1) * P(grad(v),n_phi1) + u * v, definedonelements=reg_Th)
        a += SymbolicBFI(form = (lam_nd * grad(u)*n_phi1) * (grad(v)*n_phi1), definedonelements=reg_Th)
        
        a += SymbolicBFI(levelset_domain = lset_if, form = ( - InnerProduct(grad(u),conormal1)*jump(v, vhat)
                                                                - InnerProduct(grad(v),conormal1)*jump(u, uhat)
                                                                + beta_E/h * (u- uhat)*(v - vhat) ) *gf_reg_Fh,
                                                                element_boundary=True, definedonelements=reg_Th)
        
        a += SymbolicBFI(form = beta_F/(h*h) * (u- uhat) * (v - vhat) *gf_reg_Fh, element_boundary=True, definedonelements=reg_Th)

        a += SymbolicBFI(form = beta_F2 * (grad(u)*normal_helper- sigmahat) * (grad(v)*normal_helper - tauhat) *gf_reg_Fh, element_boundary=True, definedonelements=reg_Th)

        f_coeff = -(4*(y - 1)**3 + 12)**2*((4.2*(x - 2)**2 + 4*(x - 1)**3)**2 + (4*(y - 1)**3 + 12)**2)**(-1.0)*sin(y) + (4*(y - 1)**3 + 12)*((4.2*(x - 2)**2 + 4*(x - 1)**3)**2 + (4*(y - 1)**3 + 12)**2)**(-0.5)*(-12.0*(y - 1)**2*(4*(y - 1)**3 + 12)**2*((4.2*(x - 2)**2 + 4*(x - 1)**3)**2 + (4*(y - 1)**3 + 12)**2)**(-1.5) + 12*(y - 1)**2*((4.2*(x - 2)**2 + 4*(x - 1)**3)**2 + (4*(y - 1)**3 + 12)**2)**(-0.5) - 0.5*(4.2*(x - 2)**2 + 4*(x - 1)**3)**2*((4.2*(x - 2)**2 + 4*(x - 1)**3)**2 + (4*(y - 1)**3 + 12)**2)**(-1.5)*(16.8*x + 24*(x - 1)**2 - 33.6) + ((4.2*(x - 2)**2 + 4*(x - 1)**3)**2 + (4*(y - 1)**3 + 12)**2)**(-0.5)*(8.4*x + 12*(x - 1)**2 - 16.8))*cos(y) + 2*sin(y)
        
        f = LinearForm(Vhg)
        f += SymbolicLFI(levelset_domain = lset_if, form = f_coeff * v, definedonelements=reg_Th)

        mesh.SetDeformation(deformation)
        a.Assemble()
        f.Assemble();
        
        gfu.vec[:] = 0.0
        if not condense:
            gfu.vec.data = a.mat.Inverse(freedofs, "sparsecholesky") * f.vec
        else:
            f.vec.data += a.harmonic_extension_trans * f.vec 
            inv = a.mat.Inverse(freedofs, "sparsecholesky")
            gfu.vec.data = inv * f.vec
            
            gfu.vec.data += a.harmonic_extension * gfu.vec
            gfu.vec.data += a.inner_solve * f.vec
        
        err_sqr_coefs = (gfu.components[0]-exact)**2
        l2error = sqrt( Integrate( levelset_domain=lset_if, cf=err_sqr_coefs, mesh=mesh, order=2*order+1) )
        print ("l2error : ", l2error)
        l2errors.append(l2error)
        
        H1_err_sqr_coefs = Norm( P (gfu.components[0].Deriv() - exact_grad, n_phi1) )**2
        h1error = sqrt( l2error**2 + Integrate( levelset_domain=lset_if, cf=H1_err_sqr_coefs, mesh=mesh, order=2*order+1) )
        print("h1error : ", h1error)
        h1errors.append(h1error)
        mesh.UnsetDeformation()
    
    assert l2error < 1e-7
    assert h1error < 1e-5

    if len(l2errors) > 1:
        eocs = [log(l2errors[i-1]/l2errors[i])/log(2) for i in range(1,len(l2errors))]
        print ("l2 eocs : ", eocs)
        assert sum(eocs)/len(eocs) > 5.5
        
        eocs = [log(h1errors[i-1]/h1errors[i])/log(2) for i in range(1,len(h1errors))]
        print ("h1 eocs : ", eocs)
        assert sum(eocs)/len(eocs) > 4.5
