"""Tests for the Algoim-based cut quadrature (cut_quadrature="algoim"/"auto").

Algoim (R. I. Saye, https://algoim.github.io/) computes quadrature rules
w.r.t. the higher order level set function itself (no piecewise linear
approximation, no isoparametric mapping).
"""
import pytest
import subprocess
import sys
from math import pi, sqrt
from ngsolve import *
from ngsolve.meshes import Make1DMesh, MakeStructured2DMesh, MakeStructured3DMesh
from xfem import *
from xfem.mlset import DomainTypeArray
ngsglobals.msg_level = 0


@pytest.fixture(autouse=True)
def default_globals():
    """Every test starts (and ends) with the default global settings."""
    ngsxfemglobals.SetDefaults()
    yield
    ngsxfemglobals.SetDefaults()


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
    assert ngsxfemglobals.cut_quadrature == CUT_QUADRATURE.AUTO
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
    lset_p1 = lset_gf(mesh, 1, x * x + y * y - 0.25)
    lset_p2 = lset_gf(mesh, 2, x * x + y * y - 0.25)
    one = CoefficientFunction(1)
    std = integrate(lset_p1, NEG, mesh, method="standard")
    alg = integrate(lset_p2, NEG, mesh, method="algoim")
    assert abs(alg - pi / 4) < 1e-4 < abs(std - pi / 4)

    ngsxfemglobals.cut_quadrature = "algoim"
    assert integrate(lset_p2, NEG, mesh) == pytest.approx(alg, abs=1e-14)
    with pytest.raises(Exception):
        integrate(lset_p2, NEG, mesh, method="standard")
    with pytest.raises(Exception):
        Integrate(one * dCut(lset_p2, NEG, order=4, cut_quadrature="standard"), mesh)
    ngsxfemglobals.cut_quadrature = "standard"
    with pytest.raises(Exception):
        integrate(lset_p2, NEG, mesh)
    assert integrate(lset_p1, NEG, mesh) == pytest.approx(std, abs=1e-14)
    assert integrate(lset_p2, NEG, mesh, method="algoim") == pytest.approx(alg, abs=1e-14)
    assert Integrate(one * dCut(lset_p2, NEG, order=4, cut_quadrature="algoim"), mesh) \
        == pytest.approx(alg, abs=1e-14)


def test_standard_requires_p1():
    """the standard rules reject higher order level sets (unless subdivlvl > 0)"""
    mesh = square_mesh(4, quads=False)
    levelset = (x - 0.03)**2 + (y + 0.02)**2 - 0.25
    lset_p2 = lset_gf(mesh, 2, levelset)
    for lset in [lset_p2, levelset]:
        with pytest.raises(Exception, match="piecewise linear"):
            integrate(lset, NEG, mesh, method="standard")
        with pytest.raises(Exception, match="piecewise linear"):
            CutInfo(mesh, lset, cut_quadrature="standard")
    # explicit subdivision approach (piecewise linear interpolation on subdivided simplices)
    val = Integrate(levelset_domain={"levelset": lset_p2, "domain_type": NEG, "subdivlvl": 2,
                                     "cut_quadrature": "standard"}, cf=1, mesh=mesh, order=2)
    assert val == pytest.approx(pi / 4, abs=3e-2)


USER_INFO_SCRIPT = """
from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh
from xfem import *
{setting}
mesh = MakeStructured2DMesh(nx=4, ny=4)
lset = GridFunction(H1(mesh, order={order}))
lset.Set((x - 0.5)**2 + (y - 0.5)**2 - 0.1)
Integrate(lset * dCut(lset, NEG), mesh)
Integrate(lset * dCut(lset, NEG), mesh)
"""


@pytest.mark.parametrize("order, setting, info", [(2, "", True), (1, "", False),
                                                  (2, "ngsxfemglobals.cut_quadrature = 'auto'", False),
                                                  (2, "ngsxfemglobals.cut_quadrature = 'algoim'", False)])
def test_user_info(order, setting, info):
    """the automatic choice of Algoim is reported (once), unless chosen explicitly"""
    script = USER_INFO_SCRIPT.format(order=order, setting=setting)
    out = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, check=True).stdout
    assert out.count("Algoim-based cut quadrature is used") == (1 if info else 0)


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
# facet (skeleton), element boundary and boundary cut integrals
# ---------------------------------------------------------------------------

def facet_integrals(mesh, lset, method, order=8):
    """measures of the cut parts of interior facets (skeleton), of element boundaries
    and of the domain boundary for the domain types NEG and IF"""
    V = L2(mesh, order=0, dgjumps=True)
    u, v = V.TnT()
    one = GridFunction(V)
    one.Set(1)
    res = {}
    for dt in [NEG, IF]:
        for name, form, kwargs in [("skel", 0.25 * (u + u.Other()) * (v + v.Other()), dict(skeleton=True)),
                                   ("eb", u * v, dict(element_boundary=True))]:
            a = BilinearForm(V, check_unused=False)
            a += form * dCut(lset, dt, order=order, cut_quadrature=method, **kwargs)
            a.Assemble()
            res[name, dt] = InnerProduct(one.vec, a.mat * one.vec)
        res["bnd", dt] = Integrate(CoefficientFunction(1) * dCut(lset, dt, vb=BND, order=order,
                                                                 cut_quadrature=method), mesh)
    return res


def exact_facet_integrals(mesh, cut):
    """cut(vertices of a facet) -> (measure of the NEG part, measure of the IF part)"""
    import numpy as np
    res = {(k, dt): 0.0 for k in ["skel", "eb", "bnd"] for dt in [NEG, IF]}
    for f in (mesh.edges if mesh.dim == 2 else mesh.faces):
        verts = [np.array(mesh[v].point) for v in mesh[f].vertices]
        nel = len(mesh[f].elements)
        for dt, val in zip([NEG, IF], cut(verts)):
            res["skel" if nel == 2 else "bnd", dt] += val
            res["eb", dt] += nel * val
    return res


@pytest.mark.parametrize("quads", [False, True])
def test_facets_2d(quads):
    """circle (P2 level set): the restriction to edges is exact"""
    import numpy as np
    r, c = 0.5, np.array([0.03, -0.02])

    def cut(verts):
        P, Q = verts
        d = Q - P
        a, b, cc = d @ d, 2 * d @ (P - c), (P - c) @ (P - c) - r * r
        disc = b * b - 4 * a * cc
        if disc <= 0:
            return 0.0, 0
        t1, t2 = (-b - sqrt(disc)) / (2 * a), (-b + sqrt(disc)) / (2 * a)
        return max(min(t2, 1) - max(t1, 0), 0) * sqrt(a), sum(0 < t < 1 for t in (t1, t2))

    mesh = square_mesh(8, quads, a=0.75)
    lset = lset_gf(mesh, 2, (x - c[0])**2 + (y - c[1])**2 - r**2)
    res, ref = facet_integrals(mesh, lset, "algoim"), exact_facet_integrals(mesh, cut)
    for key in ref:
        assert res[key] == pytest.approx(ref[key], abs=1e-12)
    assert ref["bnd", IF] == 0 and ref["skel", IF] > 0
    # boundary interface (points on the boundary of [-1,1]^2)
    mesh = square_mesh(4, quads)
    lset = lset_gf(mesh, 2, x * x + y * y - 1.2)
    assert Integrate(x * x * y * y * dCut(lset, IF, vb=BND), mesh) == pytest.approx(8 * 0.2, abs=1e-13)


@pytest.mark.parametrize("method", ["standard", "algoim"])
@pytest.mark.parametrize("hexes", [False, True])
def test_facets_3d_planar(hexes, method):
    """planar cut: facet, element boundary and boundary integrals are exact"""
    import numpy as np
    coefs = (0.31, -1.0, 0.17, 0.1037)
    phi = lambda p: coefs[0] * p[0] + coefs[1] * p[1] + coefs[2] * p[2] + coefs[3]

    def cut(verts):
        if len(verts) == 4:  # quadrilateral: order the vertices along its boundary
            cen = sum(verts) / 4
            e1 = verts[0] - cen
            e2 = np.cross(np.cross(verts[1] - verts[0], verts[2] - verts[0]), e1)
            verts = sorted(verts, key=lambda v: np.arctan2((v - cen) @ e2, (v - cen) @ e1))
        neg, ifpts = [], []
        for P, Q in zip(verts, verts[1:] + verts[:1]):
            if phi(P) < 0:
                neg.append(P)
            if (phi(P) < 0) != (phi(Q) < 0):
                X = P + phi(P) / (phi(P) - phi(Q)) * (Q - P)
                neg.append(X)
                ifpts.append(X)
        area = 0.5 * np.linalg.norm(sum(np.cross(A, B) for A, B in zip(neg, neg[1:] + neg[:1]))) \
            if len(neg) > 2 else 0.0
        return area, (np.linalg.norm(ifpts[0] - ifpts[1]) if len(ifpts) == 2 else 0.0)

    mesh = cube_mesh(3, hexes)
    lset = lset_gf(mesh, 1 if method == "standard" else 2, phi((x, y, z)))
    res, ref = facet_integrals(mesh, lset, method, order=2), exact_facet_integrals(mesh, cut)
    for key in ref:
        assert res[key] == pytest.approx(ref[key], abs=1e-12)


def test_facets_3d_sphere():
    """sphere (P2 level set) on hexahedra: the interface on (interior and boundary)
    faces converges to the exact arcs"""
    import numpy as np
    r, c, n = 1.2, np.array([0.031, -0.022, 0.013]), 8

    def arc(rho, cc, lo, hi):
        # length of the circle (center cc, radius rho) inside the rectangle [lo,hi]
        angles = [0.0, 2 * pi]
        for d in range(2):
            for bound in (lo[d], hi[d]):
                s = (bound - cc[d]) / rho
                if abs(s) <= 1:
                    a0 = np.arcsin(s) if d == 1 else np.arccos(s)
                    angles += [a % (2 * pi) for a in (a0, -a0, pi - a0, pi + a0)]
        angles.sort()
        return sum(rho * (b - a) for a, b in zip(angles, angles[1:])
                   if np.all(np.abs(cc + rho * np.array([np.cos((a + b) / 2), np.sin((a + b) / 2)])
                                    - (lo + hi) / 2) < (hi - lo) / 2))

    def cut(verts):
        verts = np.array(verts)
        d = int(np.argmin(np.ptp(verts, axis=0)))  # normal direction of the face
        other = [k for k in range(3) if k != d]
        dist = verts[0, d] - c[d]
        if abs(dist) >= r:
            return 0.0, 0.0
        return 0.0, arc(sqrt(r * r - dist * dist), c[other], verts[:, other].min(axis=0),
                        verts[:, other].max(axis=0))

    mesh = cube_mesh(n, hexes=True)
    lset = lset_gf(mesh, 2, (x - c[0])**2 + (y - c[1])**2 + (z - c[2])**2 - r**2)
    res, ref = facet_integrals(mesh, lset, "auto"), exact_facet_integrals(mesh, cut)
    for key in ["skel", "eb", "bnd"]:
        assert res[key, IF] == pytest.approx(ref[key, IF], abs=1e-5)
    assert ref["bnd", IF] > 0


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

    # with mesh deformation: P1 level set (isoparametric approach) -> standard rules,
    # higher order level set -> algoim
    one = CoefficientFunction(1)
    deform = GridFunction(VectorH1(mesh, order=2))
    deform.Set((0.01 * x * y, 0.02 * x))
    mesh.SetDeformation(deform)
    std_deformed = integrate(lset_p1, NEG, mesh, method="standard")
    alg_deformed = integrate(lset_p2, NEG, mesh, method="algoim")
    assert integrate(lset_p1, NEG, mesh, method="auto") == pytest.approx(std_deformed, abs=1e-14)
    assert integrate(lset_p2, NEG, mesh, method="auto") == pytest.approx(alg_deformed, abs=1e-14)
    assert abs(alg_deformed - alg) > 1e-6
    mesh.UnsetDeformation()
    # same with the deformation argument of dCut
    assert Integrate(one * dCut(lset_p1, NEG, order=4, deformation=deform), mesh) \
        == pytest.approx(std_deformed, abs=1e-14)
    assert Integrate(one * dCut(lset_p2, NEG, order=4, deformation=deform), mesh) \
        == pytest.approx(alg_deformed, abs=1e-14)


@pytest.mark.parametrize("quads", [False, True])
def test_nonaffine_deformation(quads):
    """P2 level set (exact circle in reference coordinates) on a mesh with the
    non-affine deformation Phi(x,y) = (x + a y^2, y + b x^2) (exactly represented in
    VectorH1 of order 2): the rules for the deformed domain converge to the exact
    area and length of Phi(disk) and Phi(circle)"""
    import numpy as np
    r, c, a, b = 0.5, (0.03, -0.02), 0.2, 0.15
    # det DPhi = 1 - 4ab xy, hence |Phi(disk)| = pi r^2 (1 - 4ab cx cy)
    area = pi * r**2 * (1 - 4 * a * b * c[0] * c[1])
    # length of Phi(circle): periodic trapezoidal rule (spectrally accurate)
    t = np.linspace(0, 2 * pi, 4000, endpoint=False)
    gx, gy = c[0] + r * np.cos(t), c[1] + r * np.sin(t)
    dgx, dgy = -r * np.sin(t), r * np.cos(t)
    length = 2 * pi * np.mean(np.hypot(dgx + 2 * a * gy * dgy, dgy + 2 * b * gx * dgx))
    one = CoefficientFunction(1)
    errs = []
    for n in [8, 16]:
        mesh = square_mesh(n, quads)
        lset = lset_gf(mesh, 2, (x - c[0])**2 + (y - c[1])**2 - r**2)
        deform = GridFunction(VectorH1(mesh, order=2))
        deform.Set((a * y * y, b * x * x))
        # default ("auto"): Algoim for the P2 level set, also with deformation
        neg = Integrate(one * dCut(lset, NEG, order=8, deformation=deform), mesh)
        itf = Integrate(one * dCut(lset, IF, order=8, deformation=deform), mesh)
        errs.append([abs(neg - area), abs(itf - length)])
        # same with mesh.SetDeformation and Integrate with levelset_domain
        mesh.SetDeformation(deform)
        assert integrate(lset, NEG, mesh, order=8) == pytest.approx(neg, abs=1e-13)
        assert integrate(lset, IF, mesh, order=8) == pytest.approx(itf, abs=1e-13)
        mesh.UnsetDeformation()
    for i, tol in enumerate([1e-10, 1e-8]):
        assert errs[-1][i] < tol
        assert errs[-1][i] < 1e-2 * errs[0][i]


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
