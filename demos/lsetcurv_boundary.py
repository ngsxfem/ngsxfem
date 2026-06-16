"""
Isoparametric geometry when the level set crosses the domain boundary.

The standard mesh adaptation assumes the zero level set stays inside the domain;
if it *crosses* the boundary, the deformation pushes the boundary mesh off the
domain and the cut geometry loses its higher order.  Passing
``boundary_tangential=True`` (optionally a boundary name / regex / list / Region
to restrict it) keeps the deformation tangential to the boundary (u.n = 0) and
recovers the O(h^(k+1)) accuracy.  The mesh needs a sound boundary
parametrisation for ``specialcf.normal`` (e.g. a netgen mesh).

Example: a circle crossing the bottom edge of [-1,1]^2 at 45 degrees -- we report
the cut-area error and the (machine-zero) total-area defect that certifies the
boundary is preserved.
"""
from math import pi, sin, cos, radians, log

from ngsolve import CoefficientFunction, Mesh, sqrt, x, y, Integrate as NgsIntegrate
from netgen.geom2d import SplineGeometry
from xfem import NEG, Integrate
from xfem.lsetcurv import LevelSetMeshAdaptation

order = 3
Rc, theta = 0.8, radians(45.0)
yc = -1.0 + Rc * cos(theta)
levelset = sqrt(x * x + (y - yc) ** 2) - Rc
area_exact = Rc * Rc * (pi - theta + 0.5 * sin(2 * theta))   # area of {phi<0}

print(f" circle crossing the bottom edge at 45 deg, order k={order}")
print(f"{'h':>8} {'cut-area err':>14} {'eoc':>6} {'tot-area defect':>16}")
prev = None
for n in (8, 16, 32, 64):
    geo = SplineGeometry()
    geo.AddRectangle((-1, -1), (1, 1), bcs=["bottom", "right", "top", "left"])
    mesh = Mesh(geo.GenerateMesh(maxh=2.0 / n))

    adap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.1,
                                  boundary_tangential=True)
    deform = adap.CalcDeformation(levelset)

    area = Integrate(levelset_domain={"levelset": adap.lset_p1,
                                      "domain_type": NEG, "order": 2 * order + 2},
                     cf=CoefficientFunction(1.0), mesh=mesh, deformation=deform)
    err = abs(area - area_exact)

    # u.n = 0 keeps the deformed domain the rectangle, so its area stays 4
    mesh.SetDeformation(deform)
    tot_defect = abs(NgsIntegrate(CoefficientFunction(1.0), mesh, order=2 * order) - 4.0)
    mesh.UnsetDeformation()

    eoc = "" if prev is None else f"{log(prev / err) / log(2):.2f}"
    prev = err
    print(f"{2.0/n:8.4f} {err:14.3e} {eoc:>6} {tot_defect:16.2e}")
