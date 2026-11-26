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


# ---------------------------------------------------------------------------
# IF+BND (co-dimension-2) tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cut", [0.3, 0.5, 0.7])
def test_integrate_if_bnd_2d_count(cut):
    """Integrate(1 * dCut(lset, IF, vb=BND)) counts boundary–interface intersections.

    With lset = x - cut the interface {lset=0} is the vertical line x=cut.
    It intersects the boundary of the unit square in exactly 2 points:
    (cut, 0) on the bottom edge and (cut, 1) on the top edge.
    The point measure of 2 discrete points is 2.
    """
    mesh = make_2d_mesh()
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - cut, lset)

    val = Integrate(CoefficientFunction(1) * dCut(lset, IF, vb=BND), mesh)
    assert abs(val - 2.0) < 1e-10, f"Expected 2.0 (2 intersection points), got {val}"


def test_integrate_if_bnd_2d_function():
    """Integrating a function over IF+BND evaluates it at the intersection points.

    lset = x - 0.5 → intersections at (0.5, 0) and (0.5, 1).
    ∫ x dσ₀ = 0.5 + 0.5 = 1.0
    ∫ y dσ₀ = 0.0 + 1.0 = 1.0
    """
    mesh = make_2d_mesh()
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - 0.5, lset)

    val_x = Integrate(x * dCut(lset, IF, vb=BND), mesh)
    val_y = Integrate(y * dCut(lset, IF, vb=BND), mesh)
    assert abs(val_x - 1.0) < 1e-10, f"∫x dσ₀ expected 1.0, got {val_x}"
    assert abs(val_y - 1.0) < 1e-10, f"∫y dσ₀ expected 1.0, got {val_y}"


def test_integrate_if_bnd_2d_definedon():
    """definedon restricts IF+BND to a single boundary edge."""
    mesh = make_2d_mesh()
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - 0.5, lset)

    # Only the bottom edge should have one intersection point at (0.5, 0)
    val_bot = Integrate(CoefficientFunction(1) * dCut(lset, IF, vb=BND,
                        definedon=mesh.Boundaries("bottom")), mesh)
    assert abs(val_bot - 1.0) < 1e-10, f"Expected 1 intersection on bottom, got {val_bot}"

    # Only the left edge (x=0) has no intersection with x=0.5
    val_left = Integrate(CoefficientFunction(1) * dCut(lset, IF, vb=BND,
                         definedon=mesh.Boundaries("left")), mesh)
    assert abs(val_left - 0.0) < 1e-10, f"Expected 0 intersections on left, got {val_left}"


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


# ---------------------------------------------------------------------------
# IF+BND in 3-D (co-dimension 2: lines on the boundary faces)
# ---------------------------------------------------------------------------

PLANE_A, PLANE_C = np.array([0.3, 0.2, 0.93]), 0.5671   # plane misses all cube vertices


def _plane_cube_polygon(a, c):
    """Ordered vertices of the polygon {a.x = c} ∩ [0,1]^3."""
    import itertools
    V = np.array(list(itertools.product([0, 1], repeat=3)), float)
    P = []
    for i, j in itertools.combinations(range(8), 2):
        if np.sum(np.abs(V[i] - V[j])) == 1:
            fi, fj = a @ V[i] - c, a @ V[j] - c
            if fi * fj < 0:
                P.append(V[i] + fi / (fi - fj) * (V[j] - V[i]))
    P = np.array(P)
    m = P.mean(axis=0)
    n = a / np.linalg.norm(a)
    e1 = np.cross(n, np.eye(3)[np.argmin(np.abs(n))]); e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return P[np.argsort(np.arctan2((P - m) @ e2, (P - m) @ e1))]


def _polygon_line_integrals(P):
    """Perimeter and the line integral of the (linear) function x of a closed polygon."""
    seg = np.roll(P, -1, axis=0) - P
    mid = 0.5 * (np.roll(P, -1, axis=0) + P)
    lens = np.linalg.norm(seg, axis=1)
    return lens.sum(), (lens * mid[:, 0]).sum()


def _plane_lset(mesh):
    lset = GridFunction(H1(mesh, order=1))
    InterpolateToP1(PLANE_A[0] * x + PLANE_A[1] * y + PLANE_A[2] * z - PLANE_C, lset)
    return lset


def test_integrate_if_bnd_3d_length():
    """3-D: dCut(lset, IF, vb=BND) is the length measure of the contact line on the boundary.

    Planar level set in the unit cube: the line is a polygon on the cube faces, its length and
    the line integral of x are known exactly.
    """
    mesh = make_3d_mesh(maxh=0.25)
    lset = _plane_lset(mesh)
    length, int_x = _polygon_line_integrals(_plane_cube_polygon(PLANE_A, PLANE_C))

    val = Integrate(CoefficientFunction(1) * dCut(lset, IF, vb=BND), mesh)
    val_x = Integrate(x * dCut(lset, IF, vb=BND), mesh)
    assert abs(val - length) < 1e-10, f"line length expected {length}, got {val}"
    assert abs(val_x - int_x) < 1e-10, f"∫x dγ expected {int_x}, got {val_x}"


def test_integrate_if_bnd_3d_deformed():
    """3-D IF+BND with a mesh deformation (affine: faces not axis-parallel, sheared metric).

    With deformation d = A x - x the physical domain is A([0,1]^3) and the line is the image of
    the reference polygon, so the exact length is sum |A (P_k+1 - P_k)|.
    """
    mesh = make_3d_mesh(maxh=0.25)
    lset = _plane_lset(mesh)
    A = np.array([[1.1, 0.2, -0.1], [0.15, 0.9, 0.25], [-0.2, 0.1, 1.05]])
    X = (x, y, z)
    deform = GridFunction(VectorH1(mesh, order=1))
    deform.Set(CoefficientFunction(tuple(sum(A[i, j] * X[j] for j in range(3)) - X[i]
                                         for i in range(3))))
    length, int_x = _polygon_line_integrals(_plane_cube_polygon(PLANE_A, PLANE_C) @ A.T)

    dl = dCut(lset, IF, vb=BND, deformation=deform)
    val = Integrate(CoefficientFunction(1) * dl, mesh)
    val_x = Integrate(x * dl, mesh)
    assert abs(val - length) < 1e-10, f"line length expected {length}, got {val}"
    assert abs(val_x - int_x) < 1e-10, f"∫x dγ expected {int_x}, got {val_x}"


def test_integrate_if_bnd_3d_hex():
    """3-D IF+BND on a hexahedral mesh: the boundary elements are quads (planar level set, Q1)."""
    from ngsolve.meshes import MakeStructured3DMesh
    mesh = MakeStructured3DMesh(hexes=True, nx=5, ny=4, nz=3)
    lset = GridFunction(H1(mesh, order=1))
    lset.Set(PLANE_A[0] * x + PLANE_A[1] * y + PLANE_A[2] * z - PLANE_C)
    length, int_x = _polygon_line_integrals(_plane_cube_polygon(PLANE_A, PLANE_C))

    val = Integrate(CoefficientFunction(1) * dCut(lset, IF, vb=BND), mesh)
    val_x = Integrate(x * dCut(lset, IF, vb=BND), mesh)
    assert abs(val - length) < 1e-10, f"line length expected {length}, got {val}"
    assert abs(val_x - int_x) < 1e-10, f"∫x dγ expected {int_x}, got {val_x}"


@pytest.mark.parametrize("hexes", [False, True])
def test_integrate_if_bnd_3d_deformed_nonaffine(hexes):
    """3-D IF+BND with a non-affine (quadratic) deformation, as for isoparametric unfitted FEM.

    lset = z - c, d = (0, b y^2, a y (1-y)). The reference line z = c is mapped to
    (x, y + b y^2, c + a y (1-y)) on the faces x = 0 and x = 1 and to straight segments of
    length 1 on y = 0 (d = 0) and y = 1 (d = const).
    """
    from ngsolve.meshes import MakeStructured3DMesh
    a, b, c = 0.3, 0.25, 0.437
    mesh = MakeStructured3DMesh(hexes=True, nx=4, ny=4, nz=4) if hexes else make_3d_mesh(maxh=0.25)
    lset = GridFunction(H1(mesh, order=1))
    lset.Set(z - c)
    deform = GridFunction(VectorH1(mesh, order=2))
    deform.Set(CoefficientFunction((0, b * y * y, a * y * (1 - y))))
    s, w = np.polynomial.legendre.leggauss(40)
    yq, wq = 0.5 * (s + 1), 0.5 * w
    length = 2 + 2 * np.sum(wq * np.sqrt((1 + 2 * b * yq)**2 + a**2 * (1 - 2 * yq)**2))

    val = Integrate(CoefficientFunction(1) * dCut(lset, IF, vb=BND, deformation=deform, order=12), mesh)
    assert abs(val - length) < 1e-10, f"line length expected {length}, got {val}"


def test_integrate_if_bnd_3d_spacetime():
    """3-D IF+BND with a space-time level set z - c - 0.1 t: the line has length 4 at all t,
    the space-time integral over the reference time interval [0,1] is 4."""
    mesh = make_3d_mesh(maxh=0.3)
    c = 0.437
    lset = GridFunction(SpaceTimeFESpace(H1(mesh, order=1), ScalarTimeFE(1)))
    SpaceTimeInterpolateToP1(z - c - 0.1 * tref, tref, lset)
    val = Integrate(CoefficientFunction(1) * dCut(lset, IF, vb=BND, time_order=2), mesh)
    assert abs(val - 4) < 1e-10, f"space-time line measure expected 4, got {val}"


def test_integrate_if_bnd_3d_mlset():
    """3-D multi-level-set IF+BND: line z = c restricted to x < 1/2 (faces x=0, y=0, y=1)."""
    from xfem.mlset import DomainTypeArray
    mesh = make_3d_mesh(maxh=0.3)
    l1 = GridFunction(H1(mesh, order=1))
    l1.Set(z - 0.437)
    l2 = GridFunction(H1(mesh, order=1))
    l2.Set(x - 0.5)
    val = Integrate(CoefficientFunction(1) * dCut((l1, l2), DomainTypeArray((IF, NEG)), vb=BND), mesh)
    assert abs(val - 2) < 1e-10, f"line length expected 2, got {val}"
