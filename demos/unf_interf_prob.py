"""
In this example we solve a scalar *unfitted* interface problem. As a
discretisation method we use a level set based geometry description and
a Cut or Extended Finite element method with a Nitsche formulation to
impose the interface conditions.

Domain:
-------
The domain is [-1.5,1.5]^2 while the interface is described by a level
set function ( unit ball in the ||·||_4 norm ). In the discretisation
the level set function is approximated with a piecewise linear
interpolation.

PDE problem:
------------
domain equations:
    - alpha_i (u_xx + u_yy) =   f in sub-domain i, i=1,2
interface conditions:
                        [u] =    0 on interface (continuity across
                                                 the interface     ),
           [-alpha · du/dn] =    0 on interface (conservation of the
                                                 (diffusive) flux. ),
                         u  =  u_D on domain boundary.

The r.h.s. term f and the Dirichlet data u_D is chosen according to a
manufactured solution which allows us to measure errors after the
computation of a discrete solution. The coefficients alpha are
domain-wise constants which are different in the two sub-domains.

Discretisation:
---------------
* Finite element space: We consider two different formulations:
  * "CutFEM" For each sub-domain we consider the space of piecewise
    polynomials restricted to the corresponding sub-domain.
  * "XFEM" we consider a standard space of piecewise polynomials combined
    with an enrichment space.

* Variational formulation: We use a Nitsche formulation which involves
  averages of the fluxes and jumps of the solution across the interface.
  For the average we use the Hansbo-choice [1] where the average is
  adjusted to the local cut configuration in order to ensure stability
  of the resulting formulation.

Implementational aspects:
-------------------------
* Geometry approximation: To approximate the implicitly described
  geometry we use the piecewise (multi-) linear interpolation of the
  level set function. For this geometry approximation (arbitrary order)
  accurate numerical integration is provided.
  For higher order approximation an additional mesh deformation is
  applied that maps the low order interface approximation to something
  more accurate. See also jupyter tutorials `intlset` and `cutfem`

* Linear systems: A (sparse) direct solver is applied to solve the
  arising linear systems.

Literature:
-----------
[1] A. Hansbo, P. Hansbo, An unfitted finite element method, based on
    Nitsche's method, for elliptic interface problems, Comp. Meth.
    Appl. Mech. Eng., 191(47):5537-5552, 2002.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.lsetcurv import LevelSetMeshAdaptation

from math import pi

# -------------------------------- PARAMETERS ---------------------------------
# Domain corners
ll, ur = (-1.5, -1.5), (1.5, 1.5)
# Mesh size
maxh = 0.2
# Finite element space order
order = 2

# Diffusion coefficients for the sub-domains (NEG/POS):
alpha = [1.0, 2.0]
# Nitsche penalty parameter
lambda_nitsche = 20

formulation = "XFEM"  # or "CUTFEM"

# ----------------------------------- MAIN ------------------------------------
# We generate the background mesh of the domain and use a simplicity
# triangulation to obtain a mesh with quadrilaterals use
# 'quad_dominated=True'
square = SplineGeometry()
square.AddRectangle(ll, ur, bc=1)
mesh = Mesh(square.GenerateMesh(maxh=maxh, quad_dominated=False))

# Manufactured solution and corresponding r.h.s. data CoefficientFunctions:
r22 = x**2 + y**2
r44 = x**4 + y**4
r66 = x**6 + y**6
r41 = sqrt(sqrt(r44))
r4m3 = 1.0 / r41**3
solution = [1 + pi / 2 - sqrt(2.0) * cos(pi / 4 * r44), pi / 2 * r41]
coef_f = [-alpha[i] * (solution[i].Diff(x).Diff(x)
                       + solution[i].Diff(y).Diff(y)) for i in range(2)]

# Level set function of the domain (phi = ||x||_4 - 1) and its interpolation:
levelset = r41 - 1.0

if order > 1:
    lsetadap = LevelSetMeshAdaptation(mesh, order=order, levelset=levelset)
else:
    lsetadap = NoDeformation(mesh, levelset)
lsetp1 = lsetadap.lset_p1

# Background FESpaces (used as CutFESpaces later-on):
Vh = H1(mesh, order=order, dirichlet=".*")

# Gathering information on cut elements:
ci = CutInfo(mesh, lsetp1)

if formulation == "XFEM":
    Vhx = XFESpace(Vh, lsetp1)
    VhG = Vh * Vhx

    gfu = GridFunction(VhG)

    (u_std, u_x), (v_std, v_x) = VhG.TnT()

    u = [u_std + op(u_x) for op in [neg, pos]]
    v = [v_std + op(v_x) for op in [neg, pos]]
    gradu = [grad(u_std) + op(u_x) for op in [neg_grad, pos_grad]]
    gradv = [grad(v_std) + op(v_x) for op in [neg_grad, pos_grad]]

    gfu_std, gfu_x = gfu.components
    gfuh = [gfu_std + op(gfu_x) for op in [neg, pos]]
else:
    # Create CutFESpace by selecting only active elements for inside
    # and outside:
    VhG = Compress(Vh, GetDofsOfElements(Vh, ci.GetElementsOfType(HASNEG))) \
        * Compress(Vh, GetDofsOfElements(Vh, ci.GetElementsOfType(HASPOS)))

    gfu = GridFunction(VhG)

    u, v = VhG.TnT()
    gradu, gradv = [[grad(w[i]) for i in [0, 1]] for w in [u, v]]

    gfuh = gfu.components

print("unknowns in background FESpace : ", Vh.ndof)
print("unknowns in " + formulation + " FESpace : ", VhG.ndof)


# Coefficients / parameters:
n = 1.0 / grad(lsetp1).Norm() * grad(lsetp1)
h = specialcf.mesh_size
# The cut ratio extracted from the cutinfo-class
kappa = (CutRatioGF(ci), 1.0 - CutRatioGF(ci))
# Nitsche stabilization parameter:
stab = lambda_nitsche * (alpha[1] + alpha[0]) / h

average_flux_u = sum([- kappa[i] * alpha[i] * gradu[i] * n for i in [0, 1]])
average_flux_v = sum([- kappa[i] * alpha[i] * gradv[i] * n for i in [0, 1]])

# Integration domains for integration on negative/positive sub-domains
# and on the interface: Here, the integration is (geometrically) exact
# if the "levelset"-argument is a piecewise (multi-)linear function.
# We further provide a mesh deformation that is applied in the higher order
# case:
dx = tuple([dCut(lsetp1, dt, deformation=lsetadap.deform,
                 definedonelements=ci.GetElementsOfType(HAS(dt)))
            for dt in [NEG, POS]])
ds = dCut(lsetp1, IF, deformation=lsetadap.deform)

# Bilinear form for the unfitted Nitsche formulation:
a = BilinearForm(VhG, symmetric=True)
a += sum(alpha[i] * gradu[i] * gradv[i] * dx[i] for i in [0, 1])
a += (average_flux_u * (v[0] - v[1]) + average_flux_v * (u[0] - u[1])
      + stab * (u[0] - u[1]) * (v[0] - v[1])) * ds

# R.h.s.:
f = LinearForm(VhG)
f += sum(coef_f[i] * v[i] * dx[i] for i in [0, 1])

# Setting domain boundary conditions:
# The context manager lsetadap applies the mesh deformation
# in the higher order case:
with lsetadap:
    if formulation == "XFEM":
        gfu.components[0].Set(solution[1], BND)
    else:
        gfu.components[1].Set(solution[1], BND)

# Setting up matrix and vector
a.Assemble()
f.Assemble()

# Homogenization of boundary data and solution of linear system
f.vec.data -= a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(VhG.FreeDofs()) * f.vec

# Visualize levelset, interpolated levelset and discrete solution:
Draw(levelset, mesh, "levelset")
Draw(lsetp1, mesh, "levelset_P1")
# Note that standard netgen-gui visualization does not respect
# the discontinuities, they are smeared out. To see kinks or jumps
# more clearly increase the subdivision option in the visualization.
if formulation == "CUTFEM":
    DrawDC(lsetp1, gfu.components[0], gfu.components[1],
           mesh, "u", deformation=lsetadap.deform)
    Draw(IfPos(-lsetp1,
               alpha[0]*grad(gfu.components[0]),
               alpha[1]*grad(gfu.components[1])),
         mesh, "sigma", deformation=lsetadap.deform)
else:
    DrawDC(lsetp1,
           gfu.components[0]+neg(gfu.components[1]),
           gfu.components[0]+pos(gfu.components[1]),
           mesh, "u", deformation=lsetadap.deform)
    Draw(IfPos(-lsetp1,
               alpha[0]*grad(gfu.components[0])+neg_grad(gfu.components[1]),
               alpha[1]*grad(gfu.components[0])+pos_grad(gfu.components[1])),
         mesh, "sigma", deformation=lsetadap.deform)

# Computation of L2 error:
err_sqr = sum([(gfuh[i] - solution[i])**2 * dx[i].order(2 * order)
               for i in [0, 1]])
print("L2 error : ", sqrt(Integrate(err_sqr, mesh)))
