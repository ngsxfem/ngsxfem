from ngsolve import *
from xfem import *
from netgen.geom2d import SplineGeometry
from math import pi

def test_new_integrateX_via_circle_geom(quad_dominated=True):
    square = SplineGeometry()
    square.AddRectangle([0,0],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=quad_dominated))
    r=0.6

    domains = [NEG, POS, IF]

    levelset = sqrt(x*x+y*y)-r
    referencevals = { POS : 1-pi*r*r/4, NEG : pi*r*r/4, IF : r*pi/2}

    n_ref = 8
    order = 2
    errors = dict()
    eoc = dict()

    for key in domains:
        errors[key] = []
        eoc[key] = []

    for i in range(n_ref):
        V = H1(mesh,order=1)
        lset_approx = GridFunction(V)
        InterpolateToP1(levelset,lset_approx)
    
        f = CoefficientFunction(1)
    
        for key in domains:
            integral = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=f,order=order,domain_type=key,heapsize=1000000, use_saye = True)
            print("Result of Integration Reflevel ",i,", Key ",key," : ", integral)
            errors[key].append(abs(integral - referencevals[key]))

        if i < n_ref - 1:
            mesh.Refine()
        
    for key in domains:
        eoc[key] = [log(errors[key][i+1]/errors[key][i])/log(0.5) for i in range(n_ref-1)]

    print("L2-errors:", errors)
    print("experimental order of convergence (L2):", eoc)

    for key in domains:
        for i in range(2, len(eoc[key])):            
            assert eoc[key][i] > 1.8
