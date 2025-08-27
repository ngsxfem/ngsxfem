"""
Heat equation (fitted) with nonhomogeneous Dirichlet b.c. solved with a
DG-in-time space-time finite element discretisation

Domain:
-------
The domain of this examples is [0,1]x[0,1]x[0,1] (2D+time interval).

PDE problem:
------------
  u_t - (u_xx + u_yy) = f in  Omega=[0,1]^2
                    u = sin(pi * t) * cos(pi * x) * sin(pi * y) on dOmega
The r.h.s. term f is chosen according to a manufactured solution.

Discretisation:
---------------
* Discontinuous-in-time Galerkin space-time finite elements
* Weak imposition of boundary conditions using Nitsche's method

Implementational aspects:
-------------------------------------------------------
* A (sparse) direct solver is applied to solve the arising linear systems.

References:
-----------
All concepts that are used here are explained in the jupyter-tutorial
`spacetime.ipynb` (for the more involved case of moving domains).
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import unit_square
from ngsolve import *
from xfem import *
from math import pi
ngsglobals.msg_level = 1

# -------------------------------- PARAMETERS ---------------------------------
# Space finite element order
order = 2
# Time finite element order
k_t = order
# Final simulation time
tend = 0.5
# Time step
delta_t = 1 / 8
time_order = 2*k_t

# ----------------------------------- MAIN ------------------------------------
mesh = Mesh(unit_square.GenerateMesh(maxh=0.1, quad_dominated=False))
n = specialcf.normal(mesh.dim)
h = specialcf.mesh_size

V = H1(mesh, order=order, dirichlet=[], dgjumps=False)
# V = H1(mesh, order=order, dirichlet=".*")
tfe = ScalarTimeFE(k_t)
st_fes = tfe * V

# Fitted heat equation example
tnew = 0
told = Parameter(0)
t = told + delta_t * tref

u_exact = sin(pi * t) * cos(pi * x) * sin(pi * y)

coeff_f = u_exact.Diff(t) - (u_exact.Diff(x).Diff(x) + u_exact.Diff(y).Diff(y))
coeff_f = coeff_f.Compile()

gfu = GridFunction(st_fes)
u_last = CreateTimeRestrictedGF(gfu, 1)
u, v = st_fes.TnT()

dxt = delta_t * dxtref(mesh, time_order=time_order)
dxold = dmesh(mesh, tref=0)
dxnew = dmesh(mesh, tref=1)
# for Nitsche terms:
dst = delta_t * dxtref(mesh, time_order=time_order, skeleton=True, vb=BND)


def dt(u):
    return 1.0 / delta_t * dtref(u)


a = BilinearForm(st_fes, symmetric=False)
a += grad(u) * grad(v) * dxt
a += u * v * dxold
a += dt(u) * v * dxt

lam = 100
# Nitsche terms
a += (-1) * grad(u) * n * v * dst
a += (-1) * grad(v) * n * u * dst
a += lam * (1/h) * u * v * dst
a.Assemble()

f = LinearForm(st_fes)
f += coeff_f * v * dxt
f += u_last * v * dxold
# Nitsche terms
f += lam * (1/h) * u_exact * v * dst
f += -n * grad(v) * u_exact * dst
f.Assemble()

Draw(u_last, mesh, "u")

while tend - told.Get() > delta_t / 2:
    f.Assemble()
    gfu.vec.data = a.mat.Inverse(st_fes.FreeDofs(), "umfpack") * f.vec
    RestrictGFInTime(spacetime_gf=gfu, reference_time=1.0, space_gf=u_last)
    l2error = sqrt(Integrate((u_exact - gfu)**2 * dxnew, mesh))
    Redraw()
    told.Set(told.Get() + delta_t)
    print("\rt = {0:12.9f}, L2 error = {1:12.9e}".format(told.Get(), l2error))
