"""
Heat equation (fitted) with Dirichlet b.c. solved with a DG-in-time
space-time finite element discretisation

Domain:
-------
The domain of this examples is [0,1]x[0,1]x[0,1] (2D+time interval).

PDE problem:
------------
  u_t - (u_xx + u_yy) = f in  Omega=[0,1]^2
                    u = 0 on dOmega
The r.h.s. term f is chosen according to a manufactured solution.

Discretisation:
---------------
* Discontinuous-in-time Galerkin space-time finite elements

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
order = 1
# Time finite element order
k_t = 1
# Final simulation time
tend = 1.0
# Time step
delta_t = 1 / 32


# ----------------------------------- MAIN ------------------------------------
mesh = Mesh(unit_square.GenerateMesh(maxh=0.05, quad_dominated=False))

V = H1(mesh, order=order, dirichlet=".*")
tfe = ScalarTimeFE(k_t)
st_fes = tfe * V

# Fitted heat equation example
tnew = 0
told = Parameter(0)
t = told + delta_t * tref

u_exact = sin(pi * t) * sin(pi * x)**2 * sin(pi * y)**2

coeff_f = u_exact.Diff(t) - (u_exact.Diff(x).Diff(x) + u_exact.Diff(y).Diff(y))
coeff_f = coeff_f.Compile()

gfu = GridFunction(st_fes)
u_last = CreateTimeRestrictedGF(gfu, 1)
u, v = st_fes.TnT()

dxt = delta_t * dxtref(mesh, time_order=2)
dxold = dmesh(mesh, tref=0)
dxnew = dmesh(mesh, tref=1)


def dt(u):
    return 1.0 / delta_t * dtref(u)


a = BilinearForm(st_fes, symmetric=False)
a += grad(u) * grad(v) * dxt
a += u * v * dxold
a += dt(u) * v * dxt
a.Assemble()

f = LinearForm(st_fes)
f += coeff_f * v * dxt
f += u_last * v * dxold

u_last.Set(fix_tref(u_exact, 0))
Draw(u_last, mesh, "u")

while tend - told.Get() > delta_t / 2:
    f.Assemble()
    gfu.vec.data = a.mat.Inverse(st_fes.FreeDofs(), "") * f.vec
    RestrictGFInTime(spacetime_gf=gfu, reference_time=1.0, space_gf=u_last)
    l2error = sqrt(Integrate((u_exact - gfu)**2 * dxnew, mesh))
    Redraw()
    told.Set(told.Get() + delta_t)
    print("\rt = {0:12.9f}, L2 error = {1:12.9e}".format(told.Get(), l2error))
