"""Tests for boundary cut integration: dCut(lset, NEG/POS/IF, vb=BND).

Issue #8: integrate over the part of the boundary defined by a level set.
"""
import pytest
import numpy as np
from ngsolve import *
from xfem import *
ngsglobals.msg_level = 0


# ---------------------------------------------------------------------------
# 2-D helpers
# ---------------------------------------------------------------------------

def make_2d_mesh(maxh=0.15):
    from netgen.geom2d import unit_square
    return Mesh(unit_square.GenerateMesh(maxh=maxh))


def make_3d_mesh(maxh=0.3):
    from netgen.csg import unit_cube
    return Mesh(unit_cube.GenerateMesh(maxh=maxh))


# ---------------------------------------------------------------------------
# Scalar Integrate tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cut", [0.3, 0.5, 0.7])
def test_integrate_bnd_2d(cut):
    """Integrate(1 * dCut(lset, NEG, vb=BND)) over the unit-square boundary.

    With lset = x - cut the NEG region is the left part (x < cut).
    The boundary of the unit square has four edges:
      bottom (y=0, x in [0,1]),  top (y=1, x in [0,1])
      left   (x=0, y in [0,1]),  right (x=1, y in [0,1])
    NEG portion lengths:
      bottom: cut,  top: cut,  left: 1,  right: 0  →  total = 2*cut + 1
    """
    mesh = make_2d_mesh()
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - cut, lset)

    val = Integrate(CoefficientFunction(1) * dCut(lset, NEG, vb=BND), mesh)
    expected = 2 * cut + 1.0
    assert abs(val - expected) < 1e-10, f"Expected {expected}, got {val}"


@pytest.mark.parametrize("cut", [0.3, 0.7])
def test_integrate_bnd_2d_definedon(cut):
    """dCut with definedon=mesh.Boundaries restricts to a single edge."""
    mesh = make_2d_mesh()
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - cut, lset)

    val = Integrate(CoefficientFunction(1) * dCut(lset, NEG, vb=BND,
                    definedon=mesh.Boundaries("bottom")), mesh)
    assert abs(val - cut) < 1e-10, f"Expected {cut}, got {val}"


def test_integrate_bnd_2d_pos():
    """Integrate over POS part of boundary."""
    cut = 0.4
    mesh = make_2d_mesh()
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - cut, lset)

    neg = Integrate(CoefficientFunction(1) * dCut(lset, NEG, vb=BND), mesh)
    pos = Integrate(CoefficientFunction(1) * dCut(lset, POS, vb=BND), mesh)
    total = Integrate(CoefficientFunction(1), mesh, BND)
    assert abs(neg + pos - total) < 1e-10, "NEG + POS should equal full boundary"


def test_integrate_bnd_3d(cut=0.3):
    """3-D: integrate over BND faces of unit cube where x < cut."""
    mesh = make_3d_mesh()
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - cut, lset)

    val = Integrate(CoefficientFunction(1) * dCut(lset, NEG, vb=BND), mesh)
    # NEG faces: left (x=0, area 1), bottom (z=0) * cut portion, top, front, back
    # left face (x=0): area = 1  (fully NEG)
    # right face (x=1): area = 0  (fully POS)
    # four side faces (y=0, y=1, z=0, z=1): each has NEG area = cut * 1 = cut → 4*cut
    expected = 1.0 + 4.0 * cut
    assert abs(val - expected) < 1e-10, f"Expected {expected}, got {val}"


# ---------------------------------------------------------------------------
# BFI assembly tests
# ---------------------------------------------------------------------------

def test_bfi_assembly_2d():
    """Bilinear form on boundary cut domain equals NGSolve's ds on the same region."""
    cut = 0.4
    mesh = make_2d_mesh()
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - cut, lset)

    V = H1(mesh, order=1)
    u, v = V.TnT()

    gfu = GridFunction(V)
    gfu.Set(sin(x) * cos(y))

    # Pure-NEG boundary (left edge) using dCut
    a_cut = BilinearForm(V, check_unused=False)
    a_cut += u * v * dCut(lset, NEG, vb=BND, definedon=mesh.Boundaries("left"))
    a_cut.Assemble()

    # Reference: left edge is fully NEG, so dCut == ds("left")
    a_ref = BilinearForm(V, check_unused=False)
    a_ref += u * v * ds("left")
    a_ref.Assemble()

    w_cut = a_cut.mat * gfu.vec
    w_ref = a_ref.mat * gfu.vec
    diff = Norm(w_cut - w_ref)
    assert diff < 1e-10, f"Matrix-vector product difference too large: {diff}"


def test_lfi_assembly_2d():
    """Linear form on boundary cut domain."""
    cut = 0.4
    mesh = make_2d_mesh()
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - cut, lset)

    V = H1(mesh, order=1)
    v = V.TestFunction()

    f_cut = LinearForm(V)
    f_cut += v * dCut(lset, NEG, vb=BND, definedon=mesh.Boundaries("left"))
    f_cut.Assemble()

    f_ref = LinearForm(V)
    f_ref += v * ds("left")
    f_ref.Assemble()

    diff = Norm(f_cut.vec - f_ref.vec)
    assert diff < 1e-10, f"Vector difference too large: {diff}"


def test_bfi_assembly_2d_partial():
    """On bottom edge with cut at 0.4: dCut(NEG) captures only the left portion."""
    cut = 0.4
    mesh = make_2d_mesh()
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - cut, lset)

    V = H1(mesh, order=1)
    v = V.TestFunction()

    f_cut = LinearForm(V)
    f_cut += CoefficientFunction(1) * v * dCut(lset, NEG, vb=BND,
              definedon=mesh.Boundaries("bottom"))
    f_cut.Assemble()

    # Sum of assembled vector == integral of test function (=1 * v_j)  ≈  cut length
    val = sum(f_cut.vec)
    assert abs(val - cut) < 1e-2, f"Expected ~{cut}, got {val}"


def test_bfi_definedon_bitarray():
    """definedon can be supplied as a BitArray of boundary elements."""
    cut = 0.4
    mesh = make_2d_mesh()
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - cut, lset)

    ci = CutInfo(mesh, lset)
    hasneg_bnd = ci.GetElementsOfType(HASNEG, BND)

    V = H1(mesh, order=1)
    v = V.TestFunction()

    f_cut = LinearForm(V)
    f_cut += CoefficientFunction(1) * v * dCut(lset, NEG, vb=BND,
              definedonelements=hasneg_bnd)
    f_cut.Assemble()

    val = Integrate(CoefficientFunction(1) * dCut(lset, NEG, vb=BND), mesh)
    # Both should give the same total; assert consistent
    val2 = sum(f_cut.vec)
    assert abs(val2 - val) < 1e-2, f"definedonelements mismatch: {val2} vs {val}"
