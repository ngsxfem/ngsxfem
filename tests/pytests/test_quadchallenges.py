import pytest
from ngsolve import *
from ngsolve.meshes import *
from xfem import *
from xfem.lsetcurv import *
from math import pi

@pytest.mark.parametrize("quad", [True])
@pytest.mark.parametrize("order", [2,4,8])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_via_straight_cutted_quad3D(order, domain, quad):
    mesh = MakeStructured3DMesh(hexes = True, nx=2, ny=2, nz=2)    
    
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


@pytest.mark.parametrize("quad", [True])
@pytest.mark.parametrize("order", [4])
@pytest.mark.parametrize("domain", [NEG, POS])
@pytest.mark.parametrize("alpha", [0,1,2])
@pytest.mark.parametrize("dim", [x,y,z])

#integrate f(x) = dim^alpha on the geometry implied by phi(x,y,z) = 1 - 2*x - 2*y - 2*z
# for analytic solution see
# http://www.wolframalpha.com/input/?i=integrate+from+0+to+1%2F2+from+0+to+(1%2F2-x)+from+0+to+(1%2F2-x-y)+x%5Ealpha+dz+dy+dx
def test_new_integrateX_via_straight_cutted_quad3D_polynomial(order, domain, quad, alpha, dim):
    #ngsglobals.msg_level = 0
    mesh = MakeStructured3DMesh(hexes = quad, nx=5, ny=5, nz = 5)    

    levelset = 1 - 2*x- 2*y - 2*z
    val_pos = 2**(-alpha-3)/(alpha**3+6*alpha*alpha + 11*alpha+6)
    referencevals = {POS: val_pos, NEG: 1./(alpha+1) - val_pos}
    print("Val_pos: ", val_pos)
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = dim**alpha
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain, "quad_dir_policy" : OPTIMAL},
                         cf=f, mesh=mesh, order = order)
    print("Integral:", integral)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)

def test_new_integrateX_via_straight_cutted_quad3D_polynomial_zero_val_challenge(order=4, domain=POS, alpha=2, dim=x):
    ngsglobals.msg_level = 0

    mesh = MakeStructured3DMesh(hexes = True, nx=2, ny=2, nz= 2)    
        
    levelset = 1 - 2*x- 2*y - 2*z
    val_pos = 2**(-alpha-3)/(alpha**3+6*alpha*alpha + 11*alpha+6)
    referencevals = {POS: val_pos, NEG: 1./(alpha+1) - val_pos}
    print("Val_pos: ", val_pos)
    
    lset_approx = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lset_approx)
    
    f = dim**alpha
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain, "quad_dir_policy" : OPTIMAL},
                         cf=f, mesh=mesh, order = order)
    print("Integral:", integral)
    error = abs(integral - referencevals[domain])
    
    assert error < 5e-15*(order+1)*(order+1)
    
@pytest.mark.parametrize("quad", [True])
@pytest.mark.parametrize("order", [2,4,6])
@pytest.mark.parametrize("domain", [NEG, POS, IF])

def test_new_integrateX_TPMC_case_quad3D(order, domain, quad):
    mesh = MakeStructured3DMesh(hexes = True, nx=1, ny=1, nz = 1)    
    lset_approx = GridFunction(H1(mesh,order=1))
    for i,v in enumerate([-4,4,-1,-1,2,-3,5,-1]):
        lset_approx.vec[i] = v
    
    f = CoefficientFunction(1)
    
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
    print("Integral: ", integral)

    if domain == IF:
        assert integral < 10
    elif domain == NEG:
        assert abs(integral - 0.5167820912197415) < 0.75
    else:
        assert abs(integral - 0.4825797907263282) < 0.75

@pytest.mark.parametrize("quad", [True])
@pytest.mark.parametrize("order", [2,4])
@pytest.mark.parametrize("high_order", [False, True])

def test_new_integrateX_TPMC_case_quad3D2(order, quad, high_order):

    mesh = MakeStructured3DMesh(hexes = quad, nx=10, ny=10, nz = 10)    
    
    #phi = -4*(1-x)*(1-y)*(1-z) + 4*(1-x)*(1-y)*z -1*(1-x)*y*(1-z) - 1*(1-x)*y*z + 2*x*(1-y)*(1-z) -3 *x*(1-y)*z + 5 * x * y * (1-z) -1 *x *y*z
    phi = x *((7*y - 13) *z + 6) + y *(3 - 8 *z) + 8 *z - 4
    
    if high_order:
        print("Creating LevelSetMeshAdaptation class")
        lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)
        lsetp1 = lsetmeshadap.lset_p1
        deformation = lsetmeshadap.CalcDeformation(phi)
    
        mesh.SetDeformation(deformation)    
    else:
        lsetp1 = GridFunction(H1(mesh,order=1))
        lsetp1.Set(phi)
    
    f = CoefficientFunction(1)
    
    print("Doing integration")
    for domain in [POS, NEG, IF]:
        integral = Integrate(levelset_domain = { "levelset" : lsetp1, "domain_type" : domain},
                         cf=f, mesh=mesh, order = order)
        print("Integral: ", integral, " ; domain = ", domain)

        if domain == IF:
            assert abs( integral - 1.82169) < 5e-3
        elif domain == NEG:
            assert abs(integral - 0.51681) < 1e-3
        else:
            assert abs(integral - 0.48319) < 1e-3
