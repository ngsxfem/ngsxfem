"""
Imposing different boundary conditions on parts of the same boundary label
split by a level set function. This addresses the use case from ngsxfem issue #8.

Problem
-------
We solve a Poisson problem on the unit square with different boundary conditions
on the bottom edge, split by a level set:

  -Δu = f  in Ω = (0,1)²

  u = g_NEG  on Γ_NEG  (bottom edge, x < 0.5)
  u = g_POS  on Γ_POS  (bottom edge, x > 0.5)
  u = 0      on remaining boundary (left, right, top)

where the bottom edge is split at x = 0.5 by the level set lset = x - 0.5.

The boundary conditions are imposed weakly via Nitsche's method, integrating
only over the respective cut portion of the boundary with dCut(..., vb=BND).

The same technique can be applied to any boundary label on 2D and 3D meshes.

Usage
-----
  python3 bnd_cut.py           # interactive (opens webgui if available)
  python3 bnd_cut.py testmode  # checks error norm and exits
"""

import sys
from netgen.geom2d import unit_square
from ngsolve import *
from ngsolve.internal import visoptions
from xfem import *

ngsglobals.msg_level = 1
testmode = len(sys.argv) > 1 and sys.argv[1] == "testmode"

# -------------------------------- PARAMETERS ---------------------------------
maxh = 0.1
order = 2
lambda_nitsche = 20 * order * order

# ----------------------------------- MESH ------------------------------------
mesh = Mesh(unit_square.GenerateMesh(maxh=maxh))

# ----------------------- LEVEL SET AND P1 APPROXIMATION ---------------------
# Level set splitting the bottom edge: NEG side is x < 0.5
lset_expr = x - 0.5
lset_p1 = GridFunction(H1(mesh, order=1))
InterpolateToP1(lset_expr, lset_p1)

# ----------------------- MANUFACTURED SOLUTION ------------------------------
# Choose exact solution and right-hand side
exact = sin(pi * x) * sin(pi * y)
f_rhs = 2 * pi**2 * exact

# Boundary data matching the exact solution
g_neg = exact   # on bottom, x < 0.5
g_pos = exact   # on bottom, x > 0.5

# ----------------------- FE SPACE -------------------------------------------
V = H1(mesh, order=order, dirichlet="left|right|top")
u, v = V.TnT()

n_mesh = specialcf.normal(mesh.dim)
h = specialcf.mesh_size

# ----------------------- DIFFERENTIAL SYMBOLS --------------------------------
# Standard volume integral
dvol = dx

# Bottom boundary split by level set
ds_neg = dCut(lset_p1, NEG, vb=BND, definedon=mesh.Boundaries("bottom"))
ds_pos = dCut(lset_p1, POS, vb=BND, definedon=mesh.Boundaries("bottom"))

# ----------------------- BILINEAR AND LINEAR FORMS --------------------------
a = BilinearForm(V)
f = LinearForm(V)

# Bulk term
a += grad(u) * grad(v) * dvol
f += f_rhs * v * dvol

# Nitsche terms on bottom NEG part
a += -grad(u).Trace() * n_mesh * v * ds_neg
a += -u * grad(v).Trace() * n_mesh * ds_neg
a += (lambda_nitsche / h) * u * v * ds_neg
f += -g_neg * grad(v).Trace() * n_mesh * ds_neg
f += (lambda_nitsche / h) * g_neg * v * ds_neg

# Nitsche terms on bottom POS part
a += -grad(u).Trace() * n_mesh * v * ds_pos
a += -u * grad(v).Trace() * n_mesh * ds_pos
a += (lambda_nitsche / h) * u * v * ds_pos
f += -g_pos * grad(v).Trace() * n_mesh * ds_pos
f += (lambda_nitsche / h) * g_pos * v * ds_pos

a.Assemble()
f.Assemble()

# ----------------------- SOLVE ----------------------------------------------
gfu = GridFunction(V)
gfu.vec.data = a.mat.Inverse(V.FreeDofs()) * f.vec

# ----------------------- ERROR ----------------------------------------------
err_L2 = sqrt(Integrate((gfu - exact)**2 * dx, mesh))
print(f"L2 error: {err_L2:.3e}")

if testmode:
    assert err_L2 < 2e-3, f"L2 error too large: {err_L2}"
    print("testmode: PASSED")
    sys.exit(0)

# ----------------------- VISUALISE ------------------------------------------
try:
    from ngsolve.webgui import Draw
    Draw(gfu, mesh, "u")
    Draw(lset_p1, mesh, "lset")
    input("Press Enter to exit...")
except ImportError:
    pass
