"""
Tests for the boundary-crossing isoparametric mesh adaptation, i.e. the
``boundary_tangential`` flag of ``xfem.lsetcurv.LevelSetMeshAdaptation``, used
when the zero level set crosses the domain boundary.

We check that the higher-order convergence of a *domain* functional (area of
{phi<0}) is kept -- the property the plain adaptation loses near the boundary --
and that the deformed boundary stays on the domain boundary:

* rectangular (axis-aligned) boundary: u.n == 0 exactly, total area unchanged;
* curved (disk) boundary: high order recovered, boundary defect small/decreasing.
"""

from math import pi, sin, cos, radians, log, acos, sqrt as msqrt

from ngsolve import (CoefficientFunction, Mesh, sqrt, x, y,
                     Integrate as NgsIntegrate)
from netgen.geom2d import SplineGeometry
from xfem import NEG, Integrate
from xfem.lsetcurv import LevelSetMeshAdaptation

import pytest


RC = 0.8                      # interface radius (curvature)
RECT_AREA = 4.0               # area of [-1,1]^2


def _geometry(theta_deg):
    """circle of radius RC crossing the bottom edge y=-1 at angle theta."""
    th = radians(theta_deg)
    yc = -1.0 + RC * cos(th)
    ls = sqrt(x * x + (y - yc) ** 2) - RC
    area = RC * RC * (pi - th + 0.5 * sin(2 * th))     # exact area of {phi<0}
    return ls, area


def _mesh(n):
    # netgen rectangle [-1,1]^2 with named edges: specialcf.normal needs a sound
    # boundary parametrisation, which MakeStructured2DMesh does not provide (it
    # has zero-length left/right edges).  h ~ 2/n to mimic the n x n resolution.
    geo = SplineGeometry()
    geo.AddRectangle((-1, -1), (1, 1), bcs=["bottom", "right", "top", "left"])
    return Mesh(geo.GenerateMesh(maxh=2.0 / n))


@pytest.mark.parametrize("order", [1, 2, 3])
def test_boundary_crossing_convergence(order):
    levelset, area_exact = _geometry(45.0)

    cut_err, tot_err, un_defect = [], [], []
    nrefs = 5
    for ref in range(nrefs):
        mesh = _mesh(4 * 2 ** ref)
        adap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.1,
                                      boundary_tangential=True)
        deform = adap.CalcDeformation(levelset)
        lsetp1 = adap.lset_p1

        # domain functional: area of {phi<0} with the deformation applied
        area = Integrate(levelset_domain={"levelset": lsetp1,
                                          "domain_type": NEG, "order": 2 * order + 2},
                         cf=CoefficientFunction(1.0), mesh=mesh,
                         deformation=deform)
        cut_err.append(abs(area - area_exact))

        # boundary fidelity: total (deformed) domain area must stay 4
        mesh.SetDeformation(deform)
        tot = NgsIntegrate(CoefficientFunction(1.0), mesh, order=2 * order)
        mesh.UnsetDeformation()
        tot_err.append(abs(tot - RECT_AREA))

        un_defect.append(adap.CalcBoundaryNormalDefect())

    eoc = [log(a / b) / log(2) for a, b in zip(cut_err[:-1], cut_err[1:])]
    print(f"order {order}: cut_err = {cut_err}")
    print(f"order {order}: eoc     = {eoc}")
    print(f"order {order}: tot_err = {tot_err}")
    print(f"order {order}: u.n     = {un_defect}")

    # (1) the boundary stays *exactly* on the rectangle on every level
    assert max(tot_err) < 1e-10
    assert max(un_defect) < 1e-12

    # (2) higher-order convergence of the domain functional is preserved.
    #     average the eoc over the finer half of the levels.
    s = len(eoc) // 2
    avg_eoc = sum(eoc[s:]) / len(eoc[s:])
    assert avg_eoc > order + 0.5
    assert cut_err[-1] < (1e-4 if order == 1 else 1e-6)


def test_boundary_better_than_plain():
    """The plain adaptation lets the boundary leave the rectangle; the
    boundary-aware one (flag) does not."""
    levelset, _ = _geometry(45.0)
    mesh = _mesh(16)

    plain = LevelSetMeshAdaptation(mesh, order=3, threshold=0.1)
    plain.CalcDeformation(levelset)
    mesh.SetDeformation(plain.deform)
    tot_plain = abs(NgsIntegrate(CoefficientFunction(1.0), mesh, order=6)
                    - RECT_AREA)
    mesh.UnsetDeformation()

    new = LevelSetMeshAdaptation(mesh, order=3, threshold=0.1,
                                 boundary_tangential=True)
    new.CalcDeformation(levelset)
    mesh.SetDeformation(new.deform)
    tot_new = abs(NgsIntegrate(CoefficientFunction(1.0), mesh, order=6)
                  - RECT_AREA)
    mesh.UnsetDeformation()

    print("total-area error: plain =", tot_plain, " new =", tot_new)
    assert tot_plain > 1e-6           # plain corrupts the domain boundary
    assert tot_new < 1e-10            # new keeps it exact


def _disk_mesh(maxh, curve_order, R=1.0):
    geo = SplineGeometry()
    geo.AddCircle((0, 0), R, bc="circle")
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    mesh.Curve(curve_order)
    return mesh


@pytest.mark.parametrize("order", [2, 3])
def test_curved_boundary_crossing(order):
    """Interface circle crossing a curved (disk) domain boundary: the flag
    recovers higher-order cut-area convergence (lens) that the plain adaptation
    loses, and keeps the boundary defect small and decreasing."""
    cx, r1, r2 = 0.7, 0.5, 1.0
    levelset = sqrt((x - cx) ** 2 + y * y) - r1
    d = cx
    lens = (r1 * r1 * acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1))
            + r2 * r2 * acos((d * d + r2 * r2 - r1 * r1) / (2 * d * r2))
            - 0.5 * msqrt((-d + r1 + r2) * (d + r1 - r2)
                          * (d - r1 + r2) * (d + r1 + r2)))

    hs = [0.4, 0.2, 0.1, 0.05]
    cut_new, defect_new = [], []
    for h in hs:
        mesh = _disk_mesh(h, order)
        area0 = NgsIntegrate(CoefficientFunction(1.0), mesh, order=2 * order)
        adap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.1,
                                      boundary_tangential=True)
        deform = adap.CalcDeformation(levelset)
        area = Integrate(levelset_domain={"levelset": adap.lset_p1,
                                          "domain_type": NEG,
                                          "order": 2 * order + 2},
                         cf=CoefficientFunction(1.0), mesh=mesh,
                         deformation=deform)
        cut_new.append(abs(area - lens))
        mesh.SetDeformation(deform)
        tot = NgsIntegrate(CoefficientFunction(1.0), mesh, order=2 * order)
        mesh.UnsetDeformation()
        defect_new.append(abs(tot - area0))

    eoc = [log(a / b) / log(2) for a, b in zip(cut_new[:-1], cut_new[1:])]
    print(f"curved order {order}: cut_new = {cut_new}")
    print(f"curved order {order}: eoc     = {eoc}")
    print(f"curved order {order}: defect  = {defect_new}")

    # higher-order cut-area convergence is recovered (average of finer eocs)
    s = len(eoc) // 2
    assert sum(eoc[s:]) / len(eoc[s:]) > order + 0.0
    # boundary defect stays small -- far below the plain adaptation's ~1e-3
    # (on a curved boundary it is higher order, not necessarily machine zero,
    # since a discrete facet normal only approximates the curved normal)
    assert max(defect_new) < 1e-4


def _box_mesh(n):
    from ngsolve.meshes import MakeStructured3DMesh
    return MakeStructured3DMesh(hexes=False, nx=n, ny=n, nz=n,
                                mapping=lambda a, b, c: (2 * a - 1, 2 * b - 1,
                                                         2 * c - 1))


@pytest.mark.parametrize("order", [2, 3])
def test_3d_box_boundary_crossing(order):
    """3D: a sphere crossing the bottom face z=-1 of the box [-1,1]^3.  The flag
    keeps the deformed boundary on the box (u.n = 0) and recovers higher-order
    convergence of the cut volume.  The crossing stays well inside the face (the
    smooth-boundary assumption); the sharp box edges are not touched by the cut
    elements (ns starts at 6)."""
    from ngsolve import z
    Rc, th = 0.8, radians(45.0)
    zc = -1.0 + Rc * cos(th)
    levelset = sqrt(x * x + y * y + (z - zc) ** 2) - Rc
    a = Rc * (1.0 - cos(th))                        # cap height below z=-1
    v_exact = 4.0 / 3.0 * pi * Rc ** 3 - pi / 3.0 * a * a * (3 * Rc - a)
    BOX_VOL = 8.0

    ns = [6, 9, 13]
    hs = [2.0 / n for n in ns]
    cut_err, tot_err, un = [], [], []
    for n in ns:
        mesh = _box_mesh(n)
        adap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.1,
                                      boundary_tangential=True)
        deform = adap.CalcDeformation(levelset)
        vol = Integrate(levelset_domain={"levelset": adap.lset_p1,
                                         "domain_type": NEG,
                                         "order": 2 * order + 2},
                        cf=CoefficientFunction(1.0), mesh=mesh,
                        deformation=deform)
        cut_err.append(abs(vol - v_exact))
        mesh.SetDeformation(deform)
        tot = NgsIntegrate(CoefficientFunction(1.0), mesh, order=2 * order)
        mesh.UnsetDeformation()
        tot_err.append(abs(tot - BOX_VOL))
        un.append(adap.CalcBoundaryNormalDefect())

    eoc = [log(a_ / b_) / log(h0 / h1)
           for a_, b_, h0, h1 in zip(cut_err[:-1], cut_err[1:], hs[:-1], hs[1:])]
    print(f"3d order {order}: cut_err = {cut_err}")
    print(f"3d order {order}: eoc     = {eoc}")
    print(f"3d order {order}: tot_err = {tot_err}, u.n = {un}")

    # the deformed boundary stays on the box: u.n is machine zero.  (The total
    # box volume is preserved to higher order -- in 3D the determinant's
    # null-Lagrangian term is not machine zero for order 3, only small/decaying.)
    assert max(un) < 1e-12
    assert max(tot_err) < 1e-4
    # higher-order convergence of the cut volume
    assert sum(eoc) / len(eoc) > order + 0.5


def test_boundary_selection():
    """The treated boundaries can be selected by name/regex (default '.*' = all),
    a list, or a boundary Region.  The interface crosses the *bottom* edge:
    selecting 'bottom' (or '.*') keeps the boundary; selecting only 'left' leaves
    the (uncrossed there) bottom to the plain behaviour, so the domain drifts."""
    levelset, _ = _geometry(45.0)  # crosses the bottom edge

    def total_defect(sel):
        mesh = _mesh(16)
        # sel may be a callable(mesh) -> Region (to test the Region argument)
        kw = {} if sel is None else {
            "boundary_tangential": sel(mesh) if callable(sel) else sel}
        adap = LevelSetMeshAdaptation(mesh, order=3, threshold=0.1, **kw)
        deform = adap.CalcDeformation(levelset)
        mesh.SetDeformation(deform)
        d = abs(NgsIntegrate(CoefficientFunction(1.0), mesh, order=6) - RECT_AREA)
        mesh.UnsetDeformation()
        return d

    off = total_defect(None)
    all_ = total_defect(".*")
    all_true = total_defect(True)
    bottom = total_defect("bottom")
    left = total_defect("left")
    combo = total_defect(["left", "bottom"])
    region = total_defect(lambda m: m.Boundaries("bottom"))   # Region argument
    region_off = total_defect(lambda m: m.Boundaries("left"))
    print(f"off={off:.1e} all={all_:.1e} True={all_true:.1e} bottom={bottom:.1e} "
          f"left={left:.1e} combo={combo:.1e} region={region:.1e}")

    assert off > 1e-6                     # plain corrupts the boundary
    assert all_ < 1e-10 and all_true < 1e-10
    assert bottom < 1e-10                 # the crossed boundary is selected
    assert combo < 1e-10                  # list including 'bottom'
    assert left > 1e-6                    # bottom not selected -> still drifts
    assert region < 1e-10                 # Region(bottom): crossed -> kept
    assert region_off > 1e-6              # Region(left): bottom drifts
