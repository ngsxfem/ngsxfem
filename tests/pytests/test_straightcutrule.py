import pytest
from ngsolve import *
from xfem import *
from netgen.geom2d import SplineGeometry
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from math import pi
import netgen.meshing as ngm

@pytest.mark.parametrize("domain", [NEG, POS, IF])
@pytest.mark.parametrize("alpha", [2,4,8])

def test_polynomial_ET_Segm(domain, alpha):
    order = alpha
    m = ngm.Mesh()
    m.dim = 1
    nel = 1
    
    pnums = []
    for i in range(0, nel+1):
        pnums.append (m.Add (ngm.MeshPoint (ngm.Pnt(i/nel, 0, 0))))
    
    for i in range(0,nel):
        m.Add (ngm.Element1D ([pnums[i],pnums[i+1]], index=1))
    
    m.Add (ngm.Element0D (pnums[0], index=1))
    m.Add (ngm.Element0D (pnums[nel], index=2))
    mesh = Mesh(m)
    
    x_ast = 0.78522
    levelset = x_ast-x
    referencevals = {POS:pow(x_ast, alpha+1)/(alpha+1), NEG:(1-pow(x_ast, alpha+1))/(alpha+1), IF: pow(x_ast, alpha)}
    V = H1(mesh, order=1)
    lset_approx = GridFunction(V)
    #InterpolateToP1(levelset,lset_approx)
    lset_approx.Set(levelset)
    
    f = pow(x, alpha)
    
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
    square = SplineGeometry()
    square.AddRectangle([0,0],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=quad_dominated))
    r=0.6

    levelset = sqrt(x*x+y*y)-r
    referencevals = { POS : 1-pi*r*r/4, NEG : pi*r*r/4, IF : r*pi/2}

    n_ref = 8
    errors = []

    for i in range(n_ref):
        V = H1(mesh,order=1)
        lset_approx = GridFunction(V)
        InterpolateToP1(levelset,lset_approx)
    
        f = CoefficientFunction(1)
    
        integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                             cf=f, mesh=mesh, order = order)
        print("Result of Integration Reflevel ",i,", Key ",domain," : ", integral)
        errors.append(abs(integral - referencevals[domain]))

        if i < n_ref - 1:
            mesh.Refine()
        
    eoc = [log(errors[i+1]/errors[i])/log(0.5) for i in range(n_ref-1)]

    print("L2-errors:", errors)
    print("experimental order of convergence (L2):", eoc)

    mean_eoc_array = eoc[1:]
    mean_eoc = sum(mean_eoc_array)/len(mean_eoc_array)
    assert mean_eoc > 1.75

@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_via_straight_cutted_quad2D(order, domain, quad_dominated):
    square = SplineGeometry()
    square.AddRectangle([0,0],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=quad_dominated))
    
    levelset = 1 - 2*x - 2*y
    referencevals = {NEG: 7/8, POS: 1/8, IF: 1/sqrt(2)}
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)

@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_via_straight_cutted_quad3D(order, domain, quad_dominated):
    cube = OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) ).bc(1)
    geom = CSGeometry()
    geom.Add (cube)
    ngmesh = geom.GenerateMesh(maxh=1.3, quad_dominated=quad_dominated)
    mesh = Mesh(ngmesh)
    
    levelset = 1 - 2*x - 2*y - 2*z
    referencevals = { POS : 1./48, NEG : 47./48, IF : sqrt(3)/8 }
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)

@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("order", [2])
@pytest.mark.parametrize("domain", [NEG, POS, IF])
@pytest.mark.parametrize("dim", [x,y])

def test_new_integrateX_via_straight_cutted_quad2D(order, domain, quad_dominated, dim):
    square = SplineGeometry()
    square.AddRectangle([0,0],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=quad_dominated))
    
    levelset = 1 - 3*dim
    referencevals = {NEG: 2./3, POS: 1./3, IF: 1. }
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)
