"""
This is a simple demonstration on how to use the SpaceTime VTK output
for levelset-based moving domain problems in 2D. A circle is travelling
around in a squared background domain. A piecewise linear-in-space
approximation of the level set geometry is taken as an initial
approximation. A space-time mesh deformation is then applied to recover
higher order geometrical accuracy also in space-time. In this example
only a few time steps are carried out and in each time step a VTK output
is generated which allows to visualize the space-time geometries.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.solvers import *
from ngsolve.internal import *

from xfem import *
from xfem.lset_spacetime import *

from math import pi

ngsglobals.msg_level = 1

# -------------------------------- PARAMETERS ---------------------------------
# Half side-length of background square domain
A = 1.25
# Mesh size
maxh = 0.3

# Radius of disk (the geometry)
r0 = 0.5

# Polynomial order in time
k_t = 2
# Polynomial order in space
k_s = 2
# Polynomial order in time for level set approximation
lset_order_time = 1
# Integration order in time
time_order = 2

# End time
tend = 1
# Time-step
delta_t = tend / 4

# ----------------------------------- MAIN ------------------------------------
square = SplineGeometry()
square.AddRectangle([-A, -A], [A, A])
ngmesh = square.GenerateMesh(maxh=maxh, quad_dominated=False)
mesh = Mesh(ngmesh)

# Expression for the time variable:
told = Parameter(0)
t = told + delta_t * tref


# The data:
def x0(t):
    return r0 * cos(pi * t)


def y0(t):
    return r0 * sin(pi * t)


# Level set = x - t
levelset = sqrt((x - x0(t))**2 + (y - y0(t))**2) - 0.4


# Spatial FESpace for solution
fesp1 = H1(mesh, order=1, dgjumps=True)

# Time finite element (nodal!)
tfe = ScalarTimeFE(k_t)
# Space-time finite element space
st_fes = tfe * fesp1

# Unfitted heat equation example
lset_p1 = GridFunction(st_fes)

lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space=k_s,
                                                order_time=lset_order_time,
                                                threshold=0.1,
                                                discontinuous_qn=True)

ci = CutInfo(mesh, time_order=0)

vtk_out = [levelset, lset_p1, CoefficientFunction((lset_adap_st.deform[0],
                                                   lset_adap_st.deform[1], 0)),
           BitArrayCF(ci.GetElementsOfType(IF))]
vtk_out_names = ["levelset", "lset_p1", "deform", "cutelements"]

vtk = SpaceTimeVTKOutput(ma=mesh, coefs=vtk_out, names=vtk_out_names,
                         filename="spacetime_vtk_", subdivision_x=3,
                         subdivision_t=3)

while tend - told.Get() > delta_t / 2:
    SpaceTimeInterpolateToP1(levelset, tref, lset_p1)
    dfm = lset_adap_st.CalcDeformation(levelset, tref)
    ci.Update(lset_p1, time_order=0)

    vtk.Do(t_start=told.Get(), t_end=told.Get() + delta_t)
    told.Set(told.Get() + delta_t)
