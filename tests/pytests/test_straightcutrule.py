import pytest
from ngsolve.meshes import *
from ngsolve import *
from xfem import *
from math import pi

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

@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_via_circle_geom(quad_dominated, order, domain):
    r=0.6

    levelset = sqrt(x*x+y*y)-r
    referencevals = { POS : 1-pi*r*r/4, NEG : pi*r*r/4, IF : r*pi/2}

    n_ref = 8
    errors = []

    for i in range(n_ref):
        mesh = MakeStructured2DMesh(quads = quad_dominated, nx=2**i, ny=2**i)    

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
    
@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])
@pytest.mark.parametrize("N", [1,10,30])

def test_new_integrateX_via_straight_cutted_quad2D(order, domain, quad_dominated, N):
    mesh = MakeStructured2DMesh(quads = quad_dominated, nx=N, ny=N)    
    
    levelset = 1 - 2*x - 2*y
    referencevals = {NEG: 7/8, POS: 1/8, IF: 1/sqrt(2)}
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)

@pytest.mark.parametrize("quad_dominated", [False, True])
@pytest.mark.parametrize("order", [2])
@pytest.mark.parametrize("domain", [IF, NEG, POS])
@pytest.mark.parametrize("dim", [x,y])
@pytest.mark.parametrize("eps", [1e-1, 1e-2, 5e-3, 1e-3, 0])

def test_new_integrateX_via_orth_cutted_quad2D_epsiloned(order, domain, quad_dominated, dim, eps):
    mesh = MakeStructured2DMesh(quads = quad_dominated, nx=1, ny=1)    
    
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


@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS])
@pytest.mark.parametrize("alpha", [0,1,2])
@pytest.mark.parametrize("dim", [x,y])

#integrate f(x) = dim^alpha on the geometry implied by phi(x,y,z) = 1 - 2*x - 2*y
# for analytic solution see
# http://www.wolframalpha.com/input/?i=integrate+from+0+to+1%2F2+from+0+to+(1%2F2-x)+x%5Ealpha+dy+dx
def test_new_integrateX_via_straight_cutted_quad2D_polynomial(order, domain, quad_dominated, alpha, dim):
    mesh = MakeStructured2DMesh(quads = quad_dominated, nx=1, ny=1)    
    # square = SplineGeometry()
    # square.AddRectangle([0,0],[1,1],bc=1)
    # mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=quad_dominated))
    
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

@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_via_straight_cutted_quad3D(order, domain, quad_dominated):
    mesh = MakeStructured3DMesh(hexes = quad_dominated, nx=1, ny=1, nz=1)    
    
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

@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("order", [4])
@pytest.mark.parametrize("domain", [NEG, POS])
@pytest.mark.parametrize("alpha", [0,1,2])
@pytest.mark.parametrize("dim", [x,y,z])

#integrate f(x) = dim^alpha on the geometry implied by phi(x,y,z) = 1 - 2*x - 2*y - 2*z
# for analytic solution see
# http://www.wolframalpha.com/input/?i=integrate+from+0+to+1%2F2+from+0+to+(1%2F2-x)+from+0+to+(1%2F2-x-y)+x%5Ealpha+dz+dy+dx
def test_new_integrateX_via_straight_cutted_quad3D_polynomial(order, domain, quad_dominated, alpha, dim):
    ngsglobals.msg_level = 0
    mesh = MakeStructured3DMesh(hexes = quad_dominated, nx=1, ny=1, nz=1)    
        
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

@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("order", [2])
@pytest.mark.parametrize("domain", [NEG, POS, IF])
@pytest.mark.parametrize("dim", [x,y])

def test_new_integrateX_via_orth_cutted_quad2D(order, domain, quad_dominated, dim):
    mesh = MakeStructured2DMesh(quads = quad_dominated, nx=1, ny=1)    
    
    levelset = 1 - 3*dim
    referencevals = {NEG: 2./3, POS: 1./3, IF: 1. }
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)
    

@pytest.mark.parametrize("quad_dominated", [False, True])
@pytest.mark.parametrize("order", [2])
@pytest.mark.parametrize("domain", [NEG, POS])
@pytest.mark.parametrize("dim", [x,y,z])

def test_new_integrateX_via_orth_cutted_quad3D(order, domain, quad_dominated, dim):
    mesh = MakeStructured3DMesh(hexes = quad_dominated, nx=1, ny=1, nz=1)    
    
    levelset = 1 - 2*dim
    referencevals = { POS : 1./2, NEG : 1./2, IF : 1. }
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)
