"""Tests for the Algoim-based cut quadrature (cut_quadrature="algoim"/"auto").

Algoim (R. I. Saye, https://algoim.github.io/) computes quadrature rules
w.r.t. the higher order level set function itself (no piecewise linear
approximation, no isoparametric mapping).
"""
import pytest
from contextlib import contextmanager
from math import pi, sqrt
from ngsolve import *
from ngsolve.meshes import Make1DMesh, MakeStructured2DMesh, MakeStructured3DMesh
from xfem import *
from xfem.mlset import DomainTypeArray
ngsglobals.msg_level = 0


@pytest.fixture(autouse=True)
def restore_globals():
    """Every test starts (and ends) with the default global settings."""
    cq, order = ngsxfemglobals.cut_quadrature, ngsxfemglobals.algoim_lset_order
    ngsxfemglobals.cut_quadrature = "standard"
    ngsxfemglobals.algoim_lset_order = 4
    yield
    ngsxfemglobals.cut_quadrature = cq
    ngsxfemglobals.algoim_lset_order = order


@contextmanager
def no_simd():
    """The standard rules for higher order level sets (CF path, P1 interpolation
    on simplices) currently crash with SIMD evaluation; this is independent of
    Algoim. The comparisons with these rules use non-SIMD evaluation."""
    old = ngsxfemglobals.simd_eval
    ngsxfemglobals.SwitchSIMD(False)
    try:
        yield
    finally:
        ngsxfemglobals.SwitchSIMD(old)


def square_mesh(n, quads, a=1):
    """structured mesh of [-a,a]^2"""
    return MakeStructured2DMesh(quads=quads, nx=n, ny=n,
                                mapping=lambda x, y: (2 * a * x - a, 2 * a * y - a))


def cube_mesh(n, hexes, a=1):
    """structured mesh of [-a,a]^3"""
    return MakeStructured3DMesh(hexes=hexes, nx=n, ny=n, nz=n,
                                mapping=lambda x, y, z: (2 * a * x - a, 2 * a * y - a, 2 * a * z - a))


def lset_gf(mesh, order, cf):
    lset = GridFunction(H1(mesh, order=order))
    lset.Set(cf)
    return lset


def integrate(lset, dt, mesh, cf=1, order=4, method=None, **kwargs):
    lsetdom = {"levelset": lset, "domain_type": dt}
    if method is not None:
        lsetdom["cut_quadrature"] = method
    return Integrate(levelset_domain=lsetdom, cf=cf, mesh=mesh, order=order, **kwargs)


# ---------------------------------------------------------------------------
# options
# ---------------------------------------------------------------------------

def test_options():
    assert ngsxfemglobals.cut_quadrature == CUT_QUADRATURE.STANDARD
    for name, value in [("algoim", CUT_QUADRATURE.ALGOIM), ("Auto", CUT_QUADRATURE.AUTO),
                        ("STANDARD", CUT_QUADRATURE.STANDARD)]:
        ngsxfemglobals.cut_quadrature = name
        assert ngsxfemglobals.cut_quadrature == value
        assert CUT_QUADRATURE(name) == value
    ngsxfemglobals.cut_quadrature = CUT_QUADRATURE.ALGOIM
    assert ngsxfemglobals.cut_quadrature == CUT_QUADRATURE.ALGOIM
    with pytest.raises(Exception):
        ngsxfemglobals.cut_quadrature = "unknown"
    with pytest.raises(Exception):
        integrate(None, NEG, None, method=42)


def test_global_and_local_option():
    """the per-integral option overrides the global one"""
    mesh = square_mesh(4, quads=False)
    lset = lset_gf(mesh, 2, x * x + y * y - 0.25)
    one = CoefficientFunction(1)
    with no_simd():
        std = integrate(lset, NEG, mesh)
    alg = integrate(lset, NEG, mesh, method="algoim")
    assert abs(alg - pi / 4) < 1e-4 < abs(std - pi / 4)

    ngsxfemglobals.cut_quadrature = "algoim"
    assert integrate(lset, NEG, mesh) == pytest.approx(alg, abs=1e-14)
    with no_simd():
        assert integrate(lset, NEG, mesh, method="standard") == pytest.approx(std, abs=1e-14)
        # dCut
        assert Integrate(one * dCut(lset, NEG, order=4, cut_quadrature="standard"), mesh) \
            == pytest.approx(std, abs=1e-14)
    ngsxfemglobals.cut_quadrature = "standard"
    assert Integrate(one * dCut(lset, NEG, order=4, cut_quadrature="algoim"), mesh) \
        == pytest.approx(alg, abs=1e-14)


# ---------------------------------------------------------------------------
# planar cuts: Algoim and the standard rules are exact for polynomials
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("domain", [NEG, POS, IF])
def test_planar_segm(domain):
    mesh = Make1DMesh(3)
    lset = lset_gf(mesh, 1, 0.2 - x)
    for alpha in [0, 3, 6]:
        if domain == IF:
            continue  # no interface rules on segments (dimension 0)
        std = integrate(lset, domain, mesh, cf=x**alpha, order=alpha)
        alg = integrate(lset, domain, mesh, cf=x**alpha, order=alpha, method="algoim")
        assert alg == pytest.approx(std, abs=1e-14)


@pytest.mark.parametrize("quads", [False, True])
@pytest.mark.parametrize("domain", [NEG, POS, IF])
def test_planar_2d(quads, domain):
    # reference: standard rules on simplices (exact for polynomials on planar cuts)
    mesh_ref = square_mesh(3, quads=False)
    lset_ref = lset_gf(mesh_ref, 1, 0.3 * x - y + 0.1)
    mesh = square_mesh(3, quads)
    lset = lset_gf(mesh, 1, 0.3 * x - y + 0.1)
    for cf, order in [(1, 0), (x * x * y + y**3, 3), (x**4 * y**2, 6)]:
        ref = integrate(lset_ref, domain, mesh_ref, cf=cf, order=order)
        alg = integrate(lset, domain, mesh, cf=cf, order=order, method="algoim")
        assert alg == pytest.approx(ref, abs=1e-13)


@pytest.mark.parametrize("hexes", [False, True])
@pytest.mark.parametrize("domain", [NEG, POS, IF])
def test_planar_3d(hexes, domain):
    mesh_ref = cube_mesh(2, hexes=False)
    lset_ref = lset_gf(mesh_ref, 1, 0.3 * x - y + 0.2 * z + 0.1)
    mesh = cube_mesh(2, hexes)
    lset = lset_gf(mesh, 1, 0.3 * x - y + 0.2 * z + 0.1)
    for cf, order in [(1, 0), (x * y * z + z**2, 3)]:
        ref = integrate(lset_ref, domain, mesh_ref, cf=cf, order=order)
        alg = integrate(lset, domain, mesh, cf=cf, order=order, method="algoim")
        assert alg == pytest.approx(ref, abs=1e-13)


# ---------------------------------------------------------------------------
# curved interfaces: higher order accuracy without mesh deformation
# ---------------------------------------------------------------------------

def test_quadratic_segm():
    mesh = Make1DMesh(4)
    lset = lset_gf(mesh, 2, (x - 0.3) * (x - 0.55))  # roots inside of elements
    assert integrate(lset, NEG, mesh, method="algoim") == pytest.approx(0.25, abs=1e-14)
    assert integrate(lset, POS, mesh, method="algoim") == pytest.approx(0.75, abs=1e-14)


@pytest.mark.parametrize("quads", [False, True])
def test_circle_2d(quads):
    """exact circle (P2 level set), errors of Algoim vs. standard (P1) rules"""
    r = 0.5
    ref = {NEG: pi * r**2, POS: 4 - pi * r**2, IF: 2 * pi * r}
    errs = {dt: [] for dt in ref}
    for n in [8, 16, 32]:
        mesh = square_mesh(n, quads)
        lset = lset_gf(mesh, 2, (x - 0.03)**2 + (y + 0.02)**2 - r**2)
        lset_p1 = lset_gf(mesh, 1, (x - 0.03)**2 + (y + 0.02)**2 - r**2)
        for dt in ref:
            alg = integrate(lset, dt, mesh, order=8, method="algoim")
            std = integrate(lset_p1, dt, mesh, order=8)
            errs[dt].append(abs(alg - ref[dt]))
            # standard rules: second order geometry error
            assert abs(std - ref[dt]) > 1e2 * errs[dt][-1]
    for dt in ref:
        assert errs[dt][-1] < 1e-7
        assert errs[dt][-1] < errs[dt][0]


@pytest.mark.parametrize("hexes", [False, True])
def test_sphere_3d(hexes):
    r = 0.5
    c = (0.0123, -0.0217, 0.0311)
    ref = {NEG: 4 / 3 * pi * r**3, IF: 4 * pi * r**2}
    tol = {NEG: 1e-5, IF: 1e-3}
    mesh = cube_mesh(4, hexes)
    lset = lset_gf(mesh, 2, (x - c[0])**2 + (y - c[1])**2 + (z - c[2])**2 - r**2)
    for dt in ref:
        errs = [abs(Integrate(levelset_domain={"levelset": lset, "domain_type": dt,
                                               "cut_quadrature": "algoim", "subdivlvl": lvl},
                              cf=1, mesh=mesh, order=4) - ref[dt])
                for lvl in [0, 1]]
        # subdivlvl = 1 (helps on meshes that are coarse w.r.t. the curvature)
        assert errs[1] < tol[dt]
        assert errs[1] < errs[0]


def test_cf_levelset():
    """level set CoefficientFunction: interpolation with algoim_lset_order"""
    mesh = square_mesh(4, quads=False)
    levelset = (x - 0.03)**2 + (y + 0.02)**2 - 0.25
    gf = integrate(lset_gf(mesh, 2, levelset), NEG, mesh, method="algoim")
    ngsxfemglobals.algoim_lset_order = 2
    cf = integrate(levelset, NEG, mesh, method="algoim")
    assert cf == pytest.approx(gf, abs=1e-14)
    # non-polynomial (smooth) level set function, zero level: circle with radius 0.5
    ngsxfemglobals.algoim_lset_order = 6
    cf = integrate(exp((x - 0.03)**2 + (y + 0.02)**2) - exp(0.25), IF, mesh, order=12, method="algoim")
    assert cf == pytest.approx(pi, abs=1e-4)


def test_bnd():
    """NEG/POS integrals on boundary elements: exact for a P2 level set"""
    mesh = square_mesh(4, quads=False)
    lset = lset_gf(mesh, 2, x * x + y * y - 1.2)  # cuts the boundary of [-1,1]^2
    one = CoefficientFunction(1)
    neg = Integrate(one * dCut(lset, NEG, vb=BND, cut_quadrature="algoim"), mesh)
    pos = Integrate(one * dCut(lset, POS, vb=BND, cut_quadrature="algoim"), mesh)
    assert neg == pytest.approx(8 * sqrt(0.2), abs=1e-13)
    assert pos == pytest.approx(8 - 8 * sqrt(0.2), abs=1e-13)


# ---------------------------------------------------------------------------
# automatic selection
# ---------------------------------------------------------------------------

def test_auto():
    mesh = square_mesh(4, quads=False)
    levelset = (x - 0.03)**2 + (y + 0.02)**2 - 0.25
    lset_p1 = lset_gf(mesh, 1, levelset)
    lset_p2 = lset_gf(mesh, 2, levelset)

    # P1 level set -> standard rules
    assert integrate(lset_p1, NEG, mesh, method="auto") \
        == pytest.approx(integrate(lset_p1, NEG, mesh), abs=1e-14)
    # higher order level set -> algoim
    alg = integrate(lset_p2, NEG, mesh, method="algoim")
    assert integrate(lset_p2, NEG, mesh, method="auto") == pytest.approx(alg, abs=1e-14)
    ngsxfemglobals.algoim_lset_order = 2
    assert integrate(levelset, NEG, mesh, method="auto") == pytest.approx(alg, abs=1e-14)

    # mesh deformation (isoparametric approach) -> standard rules
    one = CoefficientFunction(1)
    deform = GridFunction(VectorH1(mesh, order=2))
    deform.Set((0.01 * x * y, 0.02 * x))
    with no_simd():
        mesh.SetDeformation(deform)
        std_deformed = integrate(lset_p2, NEG, mesh)
        assert integrate(lset_p2, NEG, mesh, method="auto") == pytest.approx(std_deformed, abs=1e-14)
        assert abs(integrate(lset_p2, NEG, mesh, method="algoim") - std_deformed) > 1e-6
        mesh.UnsetDeformation()
        # same with the deformation argument of dCut
        std_deformed = Integrate(one * dCut(lset_p2, NEG, order=4, deformation=deform), mesh)
        assert Integrate(one * dCut(lset_p2, NEG, order=4, deformation=deform, cut_quadrature="auto"),
                         mesh) == pytest.approx(std_deformed, abs=1e-14)


def test_fallback_mlset():
    """multiple level sets always use the standard rules"""
    mesh = square_mesh(4, quads=False)
    lsets = tuple(lset_gf(mesh, 1, cf) for cf in [x - 0.1, y - 0.2])
    quadrant = DomainTypeArray((NEG, NEG))
    one = CoefficientFunction(1)
    std = Integrate(one * dCut(lsets, quadrant, order=2), mesh)
    ngsxfemglobals.cut_quadrature = "algoim"
    assert Integrate(one * dCut(lsets, quadrant, order=2), mesh) == pytest.approx(std, abs=1e-14)
    assert std == pytest.approx(1.1 * 1.2, abs=1e-13)


def test_fallback_spacetime():
    """space-time integrals always use the standard rules"""
    mesh = square_mesh(2, quads=False)
    fes = SpaceTimeFESpace(H1(mesh, order=1), ScalarTimeFE(1))
    lset = GridFunction(fes)
    lset.vec[:] = -1
    vals = []
    for method in ["standard", "algoim"]:
        vals.append(Integrate(levelset_domain={"levelset": lset, "domain_type": NEG,
                                               "cut_quadrature": method},
                              cf=1, mesh=mesh, order=0, time_order=0))
    assert vals[0] == pytest.approx(4, abs=1e-13)
    assert vals[1] == pytest.approx(vals[0], abs=1e-14)


# ---------------------------------------------------------------------------
# CutInfo and assembly
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("quads", [False, True])
def test_cutinfo_and_assembly(quads):
    mesh = square_mesh(8, quads)
    r = 0.5
    levelset = (x - 0.03)**2 + (y + 0.02)**2 - r**2
    lset = lset_gf(mesh, 2, levelset)

    ci = CutInfo(mesh, lset, cut_quadrature="algoim")
    hasneg = ci.GetElementsOfType(HASNEG)
    cut = ci.GetElementsOfType(IF)
    assert sum(cut) > 0
    # all elements with interface contributions are marked as cut
    if_el = Integrate(levelset_domain={"levelset": lset, "domain_type": IF, "cut_quadrature": "algoim"},
                      cf=1, mesh=mesh, order=4, element_wise=True)
    assert all((if_el[i] > 0) == cut[i] for i in range(mesh.ne))
    # cut ratios
    neg_el = Integrate(levelset_domain={"levelset": lset, "domain_type": NEG, "cut_quadrature": "algoim"},
                       cf=1, mesh=mesh, order=4, element_wise=True)
    vol_el = Integrate(1, mesh, element_wise=True)
    ratios = ci.GetCutRatios(VOL)
    assert all(abs(ratios[i] - neg_el[i] / vol_el[i]) < 1e-12 for i in range(mesh.ne))

    # assembly with dCut: mass matrix and load vector on the cut domain
    fes = Compress(H1(mesh, order=2), GetDofsOfElements(H1(mesh, order=2), hasneg))
    u, v = fes.TnT()
    dx = dCut(lset, NEG, order=8, definedonelements=hasneg, cut_quadrature="algoim")
    ds = dCut(lset, IF, order=8, definedonelements=cut, cut_quadrature="algoim")
    m = BilinearForm(u * v * dx).Assemble()
    f = LinearForm(1 * v * dx + 1 * v * ds).Assemble()
    one = GridFunction(fes)
    one.Set(1)
    assert InnerProduct(one.vec, m.mat * one.vec) == pytest.approx(pi * r**2, abs=1e-6)
    assert InnerProduct(one.vec, f.vec) == pytest.approx(pi * r**2 + 2 * pi * r, abs=1e-5)


@pytest.mark.parametrize("simd", [True, False])
def test_simd(simd):
    """SIMD and non-SIMD evaluation give the same results"""
    mesh = square_mesh(4, quads=False)
    lset = lset_gf(mesh, 3, (x - 0.03)**2 + (y + 0.02)**2 - 0.25 + 0.1 * x**3)
    ref = integrate(lset, IF, mesh, cf=x * y + 1, order=6, method="algoim")
    old = ngsxfemglobals.simd_eval
    ngsxfemglobals.SwitchSIMD(simd)
    val = Integrate(levelset_domain={"levelset": lset, "domain_type": IF, "cut_quadrature": "algoim"},
                    cf=x * y + 1, mesh=mesh, order=6)
    ngsxfemglobals.SwitchSIMD(old)
    assert val == pytest.approx(ref, abs=1e-13)
