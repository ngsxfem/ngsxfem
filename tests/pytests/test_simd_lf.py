
from xfem import ngsxfemglobals
from timeit import default_timer as timer
import pytest
from xfem import *

from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.internal import *
from xfem import *
from xfem.lsetcurv import *

import numpy as np


@pytest.mark.parametrize("maxh", [2**(-k) for k in range(8)])
@pytest.mark.parametrize("order", [2])
def test_lf(maxh, order):

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

    Vhbase = H1(mesh, order=order, dirichlet=[], dgjumps=True)
    Vh = Restrict(Vhbase, hasneg)

    u, v = Vh.TrialFunction(), Vh.TestFunction()

    # integration domains:
    dx = dCut(lsetp1, NEG, definedonelements=hasneg, deformation=deformation)

    # R.h.s. term:
    f = LinearForm(Vh)
    f += coeff_f * v * dx
    g = LinearForm(Vh)
    g += coeff_f * v * dx

    start = timer()

    f.Assemble()
    end = timer()
    vals1 = f.vec.FV().NumPy()

    t_normal = end-start
    ngsxfemglobals.SwitchSIMD(True)
    start = timer()

    g.Assemble()
    end = timer()
    vals2 = g.vec.FV().NumPy()

    t_simd = end-start
    assert(np.linalg.norm(vals1-vals2) < 1e-8)
    # We cannot guarantee, that SIMD operations are faster then non-SIMD..
  #  assert(t_normal > t_simd)
