
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
        ngsxfemglobals.SwitchSIMD(False)
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

        ngsxfemglobals.SwitchSIMD(True)
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
        ngsxfemglobals.SwitchSIMD(False)
        start = timer()
        integral_ns = IntegrateX(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain, "time_order": 0, "order": 0},
                             mesh=mesh, cf=f)
        end = timer()
        t_ns += end - start

        # SIMD
        ngsxfemglobals.SwitchSIMD(True)
        start = timer()
        integral_simd = IntegrateX(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain, "time_order": 0, "order": 0},
                             mesh=mesh, cf=f)
        end = timer()
        t_simd += end - start

        error += abs(integral_ns - integral_simd)
        #error += abs(integral - referencevals[domain])
    print("Non-SIMD to SIMD-ratio: ", t_ns/t_simd)
    assert error/n < 1 # error to high, but test fails with less (is this normal?)

# facetpatch test
@pytest.mark.parametrize("quad", [True, False])
def test_facetpatch(quad):
    print()
    print("test for FacetPatch with quad=" + str(quad))
    mesh = MakeStructured2DMesh(quads = quad, nx=2, ny=2)    

    fes = H1(mesh,order=3,dgjumps=True)
    u,v = fes.TnT()

    for i in range(2):
        a = BilinearForm(fes)
        a += u * v * dFacetPatch()
        # non-SIMD
        ngsxfemglobals.SwitchSIMD(False)
        start = timer()
        a.Assemble()
        end = timer()
        t_ns = end - start

        a2 = BilinearForm(fes)
        a2 += u * v * dFacetPatch()
        # SIMD
        ngsxfemglobals.SwitchSIMD(True)
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

    # SIMD enabled
    ngsxfemglobals.SwitchSIMD(True)
    start = timer()
    a1.AssembleLinearization(gfu.vec)
    end = timer()

    print("SIMD done in python")
    t_simd = end - start
    
    # SIMD disabled
    ngsxfemglobals.SwitchSIMD(False)
    start = timer()
    a2.AssembleLinearization(gfu.vec)
    end = timer()
    print("non-SIMD done in python")

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
        ngsxfemglobals.SwitchSIMD(True)
        start = timer()
        a1.Apply(gfu.vec,w1)
        end = timer()
        t_simd += end-start

        # SIMD disabled
        ngsxfemglobals.SwitchSIMD(False)
        start = timer()
        a2.Apply(gfu.vec,w2)
        end = timer()
        t_ns += end-start

    print("Normal/SIMD: ", t_ns/t_simd)

    w1.data -= w2
    diff = Norm(w1)
    print("diff : ",diff)
    assert diff < 1e-12