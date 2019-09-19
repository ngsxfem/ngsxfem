import pytest
from ngsolve import *
from xfem import *
from netgen.geom2d import SplineGeometry
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from math import pi

@pytest.mark.parametrize("quad", [False])
@pytest.mark.parametrize("order", [2,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_via_circle_geom(quad, order, domain):
    square = SplineGeometry()
    square.AddRectangle([0,0],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=quad))
    r=0.6

    levelset = sqrt(x*x+y*y)-r
    referencevals = { POS : 1-pi*r*r/4, NEG : pi*r*r/4, IF : r*pi/2}

    n_ref = 8
    errors = []

    for i in range(n_ref):
        f = CoefficientFunction(1)
    
        integral = Integrate(levelset_domain = { "levelset" : levelset, "domain_type" : domain},
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

@pytest.mark.parametrize("quad", [False])
@pytest.mark.parametrize("order", [2,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_via_straight_cutted_quad2D(order, domain, quad):
    square = SplineGeometry()
    square.AddRectangle([0,0],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=quad))
    
    levelset = 1 - 2*x - 2*y
    referencevals = {NEG: 7/8, POS: 1/8, IF: 1/sqrt(2)}
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : levelset, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)

@pytest.mark.parametrize("quad", [False])
@pytest.mark.parametrize("order", [2,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_via_straight_cutted_quad3D(order, domain, quad):
    cube = OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) ).bc(1)
    geom = CSGeometry()
    geom.Add (cube)
    ngmesh = geom.GenerateMesh(maxh=1.3, quad_dominated=quad)
    mesh = Mesh(ngmesh)
    
    levelset = 1 - 2*x - 2*y - 2*z
    referencevals = { POS : 1./48, NEG : 47./48, IF : sqrt(3)/8 }
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : levelset, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)
