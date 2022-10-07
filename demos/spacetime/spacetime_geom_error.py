# ------------------------------ LOAD LIBRARIES -------------------------------
import pandas as pd
from ngsolve import *
from netgen.geom2d import SplineGeometry
from xfem import *
from math import pi
from xfem.lset_spacetime import *
ngsglobals.msg_level = 1


# -------------------------------- PARAMETERS ---------------------------------
# DISCRETIZATION PARAMETERS:
def CalcGeomError(k_t, k_s, i, j):
    n_steps = 2**i
    space_refs = j
    # Polynomial order in time for level set approximation
    lset_order_time = k_t
    # Integration order in time
    time_order = 2 * k_t
    # Time stepping parameters
    tstart = 0
    tend = 0.5
    delta_t = (tend - tstart) / n_steps
    maxh = 0.5
    # Ghost-penalty parameter
    gamma = 0.05
    # Map from reference time to physical time
    told = Parameter(tstart)
    t = told + delta_t * tref

    # PROBLEM SETUP:

    # Outer domain:
    rect = SplineGeometry()
    rect.AddRectangle([-0.6, -1], [0.6, 1])

    # Level set geometry
    # Radius of disk (the geometry)
    R = 0.5
    # Position shift of the geometry in time
    rho = (1 / (pi)) * sin(2 * pi * t)
    # Convection velocity:
    w = CoefficientFunction((0, rho.Diff(t)))
    # Level set
    r = sqrt(x**2 + (y - rho)**2)
    levelset = r - R

    # ------------------------------- MAIN ------------------------------------
    ngmesh = rect.GenerateMesh(maxh=maxh, quad_dominated=False)
    for j in range(space_refs):
        ngmesh.Refine()
    mesh = Mesh(ngmesh)
    # Space time version of Levelset Mesh Adapation object. Also offers
    # integrator helper functions that involve the correct mesh deformation
    lsetadap = LevelSetMeshAdaptation_Spacetime(mesh, order_space=k_s,
                                                order_time=lset_order_time,
                                                threshold=0.5,
                                                lset_lower_bound=-3*delta_t,
                                                lset_upper_bound=3*delta_t,
                                                discontinuous_qn=True)

    def dt(u):
        return 1.0 / delta_t * dtref(u)

    maxdist = 0

    while tend - told.Get() > delta_t / 2:
        lsetadap.CalcDeformation(levelset)
        mesh.SetDeformation(lsetadap.deform)
        maxdist = max(lsetadap.CalcMaxDistance(levelset), maxdist)
        mesh.UnsetDeformation()
        told.Set(told.Get()+delta_t)
    print(maxdist)
    return maxdist


k_t = 3
k_s = 3
df = pd.DataFrame()

for i in range(6):
    for j in range(i, i+1):
        maxdist = CalcGeomError(k_t, k_s, i, j)
        df = df.append({'i': i, 'j': j, 'maxdist': maxdist}, ignore_index=True)
df['i'] = df['i'].astype(int)
df['j'] = df['j'].astype(int)
df.to_csv("geomerror.csv", index=False)
print(df)
