'''
These test-cases compare the results of SIMD-enabled and
  SIMD-disabled evaluation of e.g. .Apply(), .Assemble(), 
  .AssembleLinearized(). We cannot guarantee that SIMD-enabled
  evaluates faster than without, therefore this is not checked; 
  but can be printed with pytest -s 'filename.py'.
last edited: 08-10-2023 (DD-MM-YYYY)
'''

from xfem import ngsxfemglobals
from timeit import default_timer as timer
import pytest
from xfem import *

from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.internal import *
from xfem import *
from xfem.lsetcurv import *
import scipy.sparse as sp

from ngsolve.meshes import *
from netgen.csg import *
from math import pi
import netgen.meshing as ngm
from xfem.lset_spacetime import *

import numpy as np


@pytest.mark.parametrize("maxh", [2**(-k) for k in range(5)])
@pytest.mark.parametrize("order", [4])
def test_lf_blf(maxh, order):
    print()
    print("test for CalcElementMatrixAdd(symboliccutbfi and CalcElementVector(symboliccutlfi) with order=" + str(order) + " and maxh=" + str(maxh))
    square = SplineGeometry()
    square.AddRectangle((-1, -1), (1, 1), bc=1)
    ngmesh = square.GenerateMesh(maxh=maxh, quad_dominated=False)
    mesh = Mesh(ngmesh)

    # Manufactured exact solution for monitoring the error
    r2 = 3 / 4  # outer radius
    r1 = 1 / 4  # inner radius
    rc = (r1 + r2) / 2.0
    rr = (r2 - r1) / 2.0
    r = sqrt(x**2 + y**2)
    levelset = IfPos(r - rc, r - rc - rr, rc - r - rr)

    exact = (20 * (r2 - sqrt(x**2 + y**2)) *
             (sqrt(x**2 + y**2) - r1)).Compile()
    coeff_f = - (exact.Diff(x).Diff(x) + exact.Diff(y).Diff(y)).Compile()

    # Higher order level set approximation
    lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.1,
                                          discontinuous_qn=True)
    deformation = lsetmeshadap.CalcDeformation(levelset)
    lsetp1 = lsetmeshadap.lset_p1

    # Element, facet and dof marking w.r.t. boundary approximation with lsetp1:
    ci = CutInfo(mesh, lsetp1)
    hasneg = ci.GetElementsOfType(HASNEG)

    Vhbase = H1(mesh, order=order, dirichlet=[], dgjumps=True, dim=mesh.dim)
    Vh = Restrict(Vhbase, hasneg)

    u, v = Vh.TrialFunction(), Vh.TestFunction()

    # integration domains:
    dx = dCut(lsetp1, NEG, definedonelements=hasneg, deformation=deformation)

    # R.h.s. term:
    f = LinearForm(Vh)
    f += CF((coeff_f,sin(x))) * v * dx
    g = LinearForm(Vh)
    g += CF((coeff_f,sin(x))) * v * dx
    a = BilinearForm(Vh)
    a += (u) * (v) * dx
    b = BilinearForm(Vh)
    b += (u) * (v) * dx
    t_normal = 0
    t_blf_normal = 0
    t_simd = 0
    t_blf_simd = 0
    n=10
    for k in range(n):
        ngsxfemglobals.simd_eval = False
        ########################
        # Time Linear Form, no simd
        start = timer()

        f.Assemble()
        end = timer()
        t_normal += end-start
        #########################
        # time BLF, no simd
        start = timer()

        a.Assemble()
        end = timer()
        t_blf_normal += end-start

        ngsxfemglobals.simd_eval = True
        #########################
        # time linear form, simd

        start = timer()

        g.Assemble()
        end = timer()

        t_simd += end-start
        #########################
        # time BLF form, simd
        start = timer()

        b.Assemble()
        end = timer()
        t_blf_simd += end-start
    t_normal /= n
    t_blf_normal /= n
    t_simd /= n
    t_blf_simd /= n

    diffvec = f.vec.CreateVector()
    diffvec.data = f.vec - g.vec
    assert(diffvec.Norm() < 1e-8)

    diffvec = a.mat.AsVector().CreateVector()
    diffvec.data = a.mat.AsVector() - b.mat.AsVector()
    assert(diffvec.Norm() < 1e-8)

    # We cannot guarantee, that SIMD operations are faster then non-SIMD..
    #assert(t_normal > t_simd)
    print("Time for linear form, no simd: ", t_normal)
    print("Time for linear form, simd: ", t_simd)
    print("ratio: ", t_normal/t_simd)
    print("Time for blf, no simd: ", t_blf_normal)
    print("Time for blf, simd: ", t_blf_simd)
    print("ratio: ", t_blf_normal/t_blf_simd)


# space-time test case
@pytest.mark.parametrize("quad", [True, False])
@pytest.mark.parametrize("domain", [NEG, POS, IF])
def test_spacetime_integrateX_via_straight_cutted_quad2Dplus1D(domain, quad):
    print()
    print("spacetime test with quad=" + str(quad) + " and domain=" + str(domain))
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
    
    f = CoefficientFunction(1)
    n = 10
    t_ns = 0
    t_simd = 0
    error = 0
    for i in range(n):
        # non-SIMD
        ngsxfemglobals.simd_eval = False
        start = timer()
        integral_ns = IntegrateX(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain, "time_order": 0, "order": 0},
                             mesh=mesh, cf=f)
        end = timer()
        t_ns += end - start

        # SIMD
        ngsxfemglobals.simd_eval = True
        start = timer()
        integral_simd = IntegrateX(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain, "time_order": 0, "order": 0},
                             mesh=mesh, cf=f)
        end = timer()
        t_simd += end - start

        error += abs(integral_ns - integral_simd)
        #error += abs(integral - referencevals[domain])
    print("Non-SIMD to SIMD-ratio: ", t_ns/t_simd)
    assert error/n < 1e-8


# facetpatch test
@pytest.mark.parametrize("quad", [True, False])
def test_facetpatch(quad):
    print()
    print("test for FacetPatch with quad=" + str(quad))
    mesh = MakeStructured2DMesh(quads = quad, nx=2, ny=2)    

    fes = H1(mesh,order=3,dgjumps=True)
    u,v = fes.TnT()

    testruns = 5

    for i in range(testruns):
        a = BilinearForm(fes)
        a += u * v * dFacetPatch()
        # non-SIMD
        ngsxfemglobals.simd_eval = False
        start = timer()
        a.Assemble()
        end = timer()
        t_ns = end - start

        a2 = BilinearForm(fes)
        a2 += u * v * dFacetPatch()
        # SIMD
        ngsxfemglobals.simd_eval = True
        start = timer()
        a2.Assemble()
        end = timer()
        t_simd = end - start

    a.mat.AsVector().data -= a2.mat.AsVector().data

    error = Norm(a.mat.AsVector())
    #error += abs(integral - referencevals[domain])
    print("Non-SIMD to SIMD-ratio: ", t_ns/t_simd)
    assert error < 1e-9


# CalcLinearizedElementMatrix-test
def test_calc_linearized():
    print()
    print("Pytest for CalcLinearizedElementMatrix")
    square = SplineGeometry()
    square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=0.8))
    
    levelset = (sqrt(sqrt(x*x*x*x+y*y*y*y)) - 1.0)

    lsetp1 = GridFunction(H1(mesh))
    lsetp1.Set(levelset)
    lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
    
    ci = CutInfo(mesh,lsetp1)
    hasneg = ci.GetElementsOfType(HASNEG)

    Vh = H1(mesh, order = 2, dirichlet=[])    
    active_dofs = GetDofsOfElements(Vh,hasneg)
    Vh = Compress(Vh,active_dofs)
    
    u,v = Vh.TnT()

    gfu = GridFunction(Vh)
    gfu.Set(sin(x))
    
    a1 = BilinearForm(Vh)
    a1 += SymbolicBFI(levelset_domain = lset_neg, form = 0.5*u**2 * v)
    a2 = BilinearForm(Vh)
    a2 += SymbolicBFI(levelset_domain = lset_neg, form = 0.5*u**2 * v)

    testruns = 5
    for i in range(testruns):
        # SIMD enabled
        ngsxfemglobals.simd_eval = True
        start = timer()
        a1.AssembleLinearization(gfu.vec)
        end = timer()

        t_simd = end - start

        # SIMD disabled
        ngsxfemglobals.simd_eval = False
        start = timer()
        a2.AssembleLinearization(gfu.vec)
        end = timer()
        
    t_ns = end - start

    print("Normal to SIMD ratio: ", str(t_ns/t_simd))

    w1 = gfu.vec.CreateVector()
    w2 = gfu.vec.CreateVector()

    w1.data = a1.mat * gfu.vec
    w2.data = a2.mat * gfu.vec
    
    # print(Norm(w1))
    # print(Norm(w2))

    w1.data -= w2.data
    diff = Norm(w1)
    print("diff : ",diff)
    assert diff < 1e-12


# test for ApplyElementMatrix
def test_calc_linearized():
    print()
    print("test for CalcLinearizedElementMatrix")
    square = SplineGeometry()
    square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=0.8))
    
    levelset = (sqrt(sqrt(x*x*x*x+y*y*y*y)) - 1.0)

    lsetp1 = GridFunction(H1(mesh))
    lsetp1.Set(levelset)
    lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
    
    ci = CutInfo(mesh,lsetp1)
    hasneg = ci.GetElementsOfType(HASNEG)

    Vh = H1(mesh, order = 2, dirichlet=[])    
    active_dofs = GetDofsOfElements(Vh,hasneg)
    Vh = Compress(Vh,active_dofs)
    
    u,v = Vh.TnT()

    gfu = GridFunction(Vh)
    gfu.Set(sin(x))

    a1 = BilinearForm(Vh)
    a1 += SymbolicBFI(levelset_domain = lset_neg, form = 0.5*u**2 * v)
    a1.Assemble()
    
    a2 = BilinearForm(Vh)
    a2 += SymbolicBFI(levelset_domain = lset_neg, form = 0.5*u**2 * v)
    a2.Assemble()

    w1 = gfu.vec.CreateVector()
    w2 = gfu.vec.CreateVector()

    w1.data = a1.mat * gfu.vec
    w2.data = a2.mat * gfu.vec
    
    t_simd = 0
    t_ns = 0
    testruns = 5
    for i in range(testruns):
        # SIMD enabled
        ngsxfemglobals.simd_eval = True
        start = timer()
        a1.Apply(gfu.vec,w1)
        end = timer()
        t_simd += end-start

        # SIMD disabled
        ngsxfemglobals.simd_eval = False
        start = timer()
        a2.Apply(gfu.vec,w2)
        end = timer()
        t_ns += end-start

    print("Normal/SIMD: ", t_ns/t_simd)

    w1.data -= w2
    diff = Norm(w1)
    print("error: ",diff)
    assert diff < 1e-12


# copy from test_straightcutrule.py - this test compares the results SIMD to non-SIMD
@pytest.mark.parametrize("i", [2, 3, 4])
def test_eb_cut_integrator_2d(i):
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

    ngsxfemglobals.simd_eval = True
    l2error_simd = sqrt( Integrate( levelset_domain=lset_if, cf=err_sqr_coefs, mesh=mesh, order=2*order+1) )
        
    ngsxfemglobals.simd_eval = False
    l2error_ns = sqrt( Integrate( levelset_domain=lset_if, cf=err_sqr_coefs, mesh=mesh, order=2*order+1) )
        
    H1_err_sqr_coefs = Norm( P (gfu.components[0].Deriv() - exact_grad, n_phi1) )**2
    ngsxfemglobals.simd_eval = True
    h1error_simd = sqrt( l2error_simd**2 + Integrate( levelset_domain=lset_if, cf=H1_err_sqr_coefs, mesh=mesh, order=2*order+1) )
    ngsxfemglobals.simd_eval = False
    h1error_ns = sqrt( l2error_ns**2 + Integrate( levelset_domain=lset_if, cf=H1_err_sqr_coefs, mesh=mesh, order=2*order+1) )

    assert abs(l2error_simd - l2error_ns) < 1e-12
    assert abs(h1error_simd - h1error_ns) < 1e-12

    mesh.UnsetDeformation()

# copy from test_differential_symbol.pytest - this test compares the result SIMD to non-SIMD
@pytest.mark.parametrize('order', [1, 2, 3])
@pytest.mark.parametrize('DOM', [POS, NEG])
@pytest.mark.parametrize('skeleton', [True, False])
@pytest.mark.parametrize('element_boundary', [True, False])
def test_cut_symbols_dg(order, DOM, skeleton, element_boundary):
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    levelset = (x - 0.5)**2 + (y - 0.5)**2 - 0.33**2

    if skeleton is False and element_boundary is False:
        return None
    if order > 1:
        lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order,
                                              levelset=levelset)
    else:
        lsetmeshadap = NoDeformation(mesh, levelset)

    lset_h = lsetmeshadap.lset_p1
    deform = lsetmeshadap.deform

    ci = CutInfo(mesh, lset_h)
    els_hasdom = ci.GetElementsOfType(HAS(DOM))
    facets_dom = GetFacetsWithNeighborTypes(mesh, a=els_hasdom, b=els_hasdom)
    els_none = BitArray(mesh.ne)
    els_none[:] = False

    V = H1(mesh, order=order, dgjumps=True)
    u, v = V.TnT()

    gfu = GridFunction(V)
    gfu.Set(sin(x))

    w1, w2 = gfu.vec.CreateVector(), gfu.vec.CreateVector()

    # Bilinear form to integrate
    nF = specialcf.normal(mesh.dim)
    flux_u = -0.5 * (grad(u) + grad(u.Other())) * nF
    flux_v = -0.5 * (grad(v) + grad(v.Other())) * nF
    jump_u = u - u.Other()
    jump_v = v - v.Other()

    form = jump_u * jump_v + flux_u * jump_v + flux_v * jump_u

    # Differential Symbol Version
    a1 = RestrictedBilinearForm(V, element_restriction=els_hasdom,
                                facet_restriction=facets_dom,
                                check_unused=False)
    a1 += form * dCut(lset_h, DOM, definedonelements=facets_dom,
                      deformation=deform, skeleton=skeleton,
                      element_boundary=element_boundary)
    a2 = RestrictedBilinearForm(V, element_restriction=els_hasdom,
                                facet_restriction=facets_dom,
                                check_unused=False)
    a2 += form * dCut(lset_h, DOM, definedonelements=facets_dom,
                      deformation=deform, skeleton=skeleton,
                      element_boundary=element_boundary)
    
    ngsxfemglobals.simd_eval = True
    a1.Assemble()
    ngsxfemglobals.simd_eval = False
    a2.Assemble()
    
    w1.data = a1.mat * gfu.vec
    w2.data = a2.mat * gfu.vec

    # # SymbolicBFI version
    # lset_dom = {'levelset': lset_h, 'domain_type': DOM, 'subdivlvl': 0}
    # a2 = RestrictedBilinearForm(V, element_restriction=els_hasdom,
    #                             facet_restriction=facets_dom,
    #                             check_unused=False)
    # a2 += SymbolicBFI(levelset_domain=lset_dom, form=form,
    #                   definedonelements=facets_dom, skeleton=skeleton,
    #                   element_boundary=element_boundary)
    # mesh.SetDeformation(deform)
    # a2.Assemble()
    # mesh.UnsetDeformation()
    # w2.data = a2.mat * gfu.vec

    # Check
    w1.data -= w2
    diff = Norm(w1)
    print(f'diff : {diff}')
    assert diff < 1e-12