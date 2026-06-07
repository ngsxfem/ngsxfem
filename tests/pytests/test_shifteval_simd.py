# Compare the SIMD and the non-SIMD evaluation path of shifted_eval.
#
# shifted_eval is realized as a (Grid)function CoefficientFunction whose
# DifferentialOperator (DiffOpShiftedEval) now provides SIMD overloads.
# Evaluating the resulting CF (e.g. via GridFunction.Set) uses the SIMD
# Apply when use_simd=True and the scalar Apply otherwise. Both paths
# compute the same shifted reference points, hence the results must agree
# up to round-off.

import pytest
from ngsolve import *
from ngsolve.meshes import *
from xfem import *

ngsglobals.msg_level = 1


def _make_mesh(dimension):
    if dimension == 2:
        return MakeStructured2DMesh(quads=False, nx=8, ny=8)
    return MakeStructured3DMesh(hexes=False, nx=8, ny=8, nz=8)


def _make_deformations(mesh, dimension, fes_dfm, use_back, use_forth):
    dfm_back = GridFunction(fes_dfm)
    dfm_forth = GridFunction(fes_dfm)
    if use_back:
        if dimension == 2:
            dfm_back.Set(CoefficientFunction((0.2 * sin(5 * y), 0.2 * cos(5 * x))))
        else:
            dfm_back.Set(CoefficientFunction((0.15 * sin(5 * y), 0.15 * cos(5 * z),
                                              0.15 * sin(5 * x))))
        # keep the deformation zero on the vertices (as in test_shifteval.py)
        for i in range(dimension * mesh.nv):
            dfm_back.vec[i] = 0.0
    if use_forth:
        if dimension == 2:
            dfm_forth.Set(CoefficientFunction((0.1 * cos(4 * x), 0.1 * sin(4 * y))))
        else:
            dfm_forth.Set(CoefficientFunction((0.1 * cos(4 * x), 0.1 * sin(4 * y),
                                               0.1 * cos(4 * z))))
        for i in range(dimension * mesh.nv):
            dfm_forth.vec[i] = 0.0
    back = dfm_back if use_back else None
    forth = dfm_forth if use_forth else None
    return back, forth


@pytest.mark.parametrize("dimension", [2, 3])
@pytest.mark.parametrize("use_back,use_forth", [(True, False), (False, True), (True, True)])
def test_shifteval_simd_vs_nosimd(dimension, use_back, use_forth):
    mesh = _make_mesh(dimension)

    fes = H1(mesh, order=3)
    fes_dfm = H1(mesh, order=3, dim=dimension)

    back, forth = _make_deformations(mesh, dimension, fes_dfm, use_back, use_forth)

    gfu_old = GridFunction(fes)
    gfu_old.Set(sin(10 * y) + cos(7 * x))

    cf = shifted_eval(gfu_old, back, forth)

    # --- Apply path: GridFunction.Set with and without SIMD ---
    gf_simd = GridFunction(fes)
    gf_nosimd = GridFunction(fes)
    gf_simd.Set(cf, use_simd=True)
    gf_nosimd.Set(cf, use_simd=False)

    diff_set = sqrt(Integrate((gf_simd - gf_nosimd) ** 2, mesh, order=8))
    ref = sqrt(Integrate(gf_nosimd ** 2, mesh, order=8))
    print(f"dim={dimension} back={use_back} forth={use_forth} "
          f"|set_simd - set_nosimd|_L2 = {diff_set:.3e} (ref {ref:.3e})")
    assert diff_set < 1e-11 * max(ref, 1.0)

    # --- AddTrans path: linear form assembly with and without SIMD ---
    v = fes.TestFunction()
    lf_simd = LinearForm(fes)
    lf_simd += SymbolicLFI(cf * v, simd_evaluate=True)
    lf_simd.Assemble()
    lf_nosimd = LinearForm(fes)
    lf_nosimd += SymbolicLFI(cf * v, simd_evaluate=False)
    lf_nosimd.Assemble()

    diffvec = lf_simd.vec.CreateVector()
    diffvec.data = lf_simd.vec - lf_nosimd.vec
    diff_lf = Norm(diffvec)
    refn = Norm(lf_nosimd.vec)
    print(f"dim={dimension} back={use_back} forth={use_forth} "
          f"|lf_simd - lf_nosimd|_2 = {diff_lf:.3e} (ref {refn:.3e})")
    assert diff_lf < 1e-11 * max(refn, 1.0)
