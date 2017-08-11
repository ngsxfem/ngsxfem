import pytest
from ngsolve import *
from xfem import *
from netgen.geom2d import SplineGeometry
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from math import pi
import netgen.meshing as ngm
from make_uniform2D_grid import MakeUniform2DGrid
from make_uniform3D_grid import MakeUniform3DGrid

@pytest.mark.parametrize("quad_dominated", [True])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_via_straight_cutted_quad3D(order, domain, quad_dominated):
    mesh = MakeUniform3DGrid(quads = True, N=2, P1=(0,0,0),P2=(1,1,1))
    
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


@pytest.mark.parametrize("quad_dominated", [True])
@pytest.mark.parametrize("order", [4])
@pytest.mark.parametrize("domain", [NEG, POS])
@pytest.mark.parametrize("alpha", [0,1,2])
@pytest.mark.parametrize("dim", [x,y,z])

#integrate f(x) = dim^alpha on the geometry implied by phi(x,y,z) = 1 - 2*x - 2*y - 2*z
# for analytic solution see
# http://www.wolframalpha.com/input/?i=integrate+from+0+to+1%2F2+from+0+to+(1%2F2-x)+from+0+to+(1%2F2-x-y)+x%5Ealpha+dz+dy+dx
def test_new_integrateX_via_straight_cutted_quad3D_polynomial(order, domain, quad_dominated, alpha, dim):
    ngsglobals.msg_level = 0

    if (quad_dominated):
        mesh = MakeUniform3DGrid(quads = True, N=2, P1=(0,0,0),P2=(1,1,1))
    else:
        cube = OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) ).bc(1)
        geom = CSGeometry()
        geom.Add (cube)
        ngmesh = geom.GenerateMesh(maxh=1.3, quad_dominated=quad_dominated)
        mesh = Mesh(ngmesh)
        
    levelset = 1 - 2*x- 2*y - 2*z
    val_pos = pow(2,-alpha-3)/(pow(alpha,3)+6*alpha*alpha + 11*alpha+6)
    referencevals = {POS: val_pos, NEG: 1./(alpha+1) - val_pos}
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = pow(dim,alpha)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)
    
