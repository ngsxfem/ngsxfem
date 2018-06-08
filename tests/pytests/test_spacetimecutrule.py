import pytest
from ngsolve import *
from xfem import *
from netgen.geom2d import SplineGeometry
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from math import pi
import netgen.meshing as ngm
from make_uniform2D_grid import MakeUniform2DGrid
from make_uniform3D_grid import MakeUniform3DGrid

tref = ReferenceTimeVariable()

@pytest.mark.parametrize("quad_dominated", [True, False])
@pytest.mark.parametrize("integrands", [(tref,0.5,0,1),
                                        (tref**3,0.25,0,3),
                                        ((1-tref)**3,0.25,0,3),
                                        (x,0.5,1,0),
                                        (tref*tref*(x*x+y*y),2/9,2,2)])
def test_spacetime_integrate_no_cut(quad_dominated, integrands):
    mesh = MakeUniform2DGrid(quads = quad_dominated, N=1, P1=(0,0), P2=(1,1))

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
    mesh = MakeUniform2DGrid(quads = quad_dominated, N=1, P1=(0,0), P2=(1,1))

    tref = ReferenceTimeVariable()
    
    levelset = lambda t : 1 - 2*x - 2*t
    referencevals = { POS : 1./8, NEG : 1 - 1/8, IF : 1.0/2 }

    h1fes = H1(mesh,order=1)
    lset_approx_h1 = GridFunction(h1fes)
    tfe = ScalarTimeFE(1) 
    fes= SpaceTimeFESpace(h1fes,tfe)
    lset_approx = GridFunction(fes)

    InterpolateToP1(levelset(0),lset_approx_h1)
    lset_approx.vec[0:h1fes.ndof] = lset_approx_h1.vec
    InterpolateToP1(levelset(1),lset_approx_h1)
    lset_approx.vec[h1fes.ndof:2*h1fes.ndof] = lset_approx_h1.vec

    print(lset_approx.vec)
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = 0, time_order=0)
    print("Integral: ", integral)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15
