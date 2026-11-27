"""
Ghost penalty (dFacetPatch) with a coefficient which is evaluated in both elements of the
patch (.Other()): the scalar evaluation (ngsxfemglobals.simd_eval = False) and the fallback
from the SIMD to the scalar evaluation (a GridFunction on another mesh has no SIMD
evaluation) give the same matrix as the SIMD evaluation.
"""
import random

import pytest
from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh
from xfem import *


def energies(kappa_on_other_mesh):
    mesh = MakeStructured2DMesh(quads=False, nx=12, ny=12)
    lsetp1 = GridFunction(H1(mesh, order=1))
    InterpolateToP1(sqrt((x - 0.5)**2 + (y - 0.5)**2) - 0.31, lsetp1)
    ci = CutInfo(mesh, lsetp1)
    gpfacets = GetFacetsWithNeighborTypes(mesh, a=ci.GetElementsOfType(HASNEG),
                                          b=ci.GetElementsOfType(IF))
    kmesh = MakeStructured2DMesh(quads=False, nx=7, ny=7) if kappa_on_other_mesh else mesh
    kappa = GridFunction(H1(kmesh, order=1))
    kappa.Set(1 + x * y)

    Vh = H1(mesh, order=2, dgjumps=True)
    u, v = Vh.TnT()
    w = GridFunction(Vh)
    random.seed(1)
    for i in range(len(w.vec)):
        w.vec[i] = random.random()
    tmp = w.vec.CreateVector()

    result = []
    for simd in [True, False]:
        ngsxfemglobals.simd_eval = simd
        try:
            a = BilinearForm(Vh)
            a += 0.5 * (kappa + kappa.Other()) * (u - u.Other()) * (v - v.Other()) \
                * dFacetPatch(definedonelements=gpfacets)
            a.Assemble()
        finally:
            ngsxfemglobals.simd_eval = True
        tmp.data = a.mat * w.vec
        result.append(InnerProduct(tmp, w.vec))
    return result


@pytest.mark.parametrize("kappa_on_other_mesh", [False, True])
def test_facetpatch_coefficient_other(kappa_on_other_mesh):
    e_simd, e_scalar = energies(kappa_on_other_mesh)
    assert e_simd > 0
    assert abs(e_simd - e_scalar) < 1e-12 * e_simd
