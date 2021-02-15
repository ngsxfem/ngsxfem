"""
In this example we solve a scalar *unfitted* interface problem. As a
discretisation method we use a level set based geometry description and
a Cut Finite element method with a Nitsche formulation to impose the
interface conditions.

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
* Finite element space: For each sub-domain we consider the space of
  piecewise linear restricted to the corresponding sub-domain.

* Variational formulation: We use a Nitsche formulation which involves
  averages of the fluxes and jumps of the solution across the interface.
  For the average we use the Hansbo-choice [1] where the average is
  adjusted to the local cut configuration in order to ensure stability
  of the resulting formulation.

Implementational aspects:
-------------------------
* Finite element space: The finite element space is a simple composition
  of two scalar piecewise linear finite element spaces w.r.t. the
  background mesh. However, in the variational formulation only these
  degrees of freedom that have some part in the corresponding sub-domain
  are used and appear in the equations. Correspondingly, degrees of
  freedom in the active part of the mesh and unused degrees of freedoms
  are marked differently to obtain a consistent treatment and a regular
  stiffness matrix.

* Geometry approximation: To approximate the implicitly described
  geometry we use the piecewise (multi-) linear interpolation of the
  level set function. For this geometry approximation (arbitrary order)
  accurate numerical integration is provided.

* Linear systems: A (sparse) direct solver is applied to solve the
  arising linear systems.

Extensions:
-----------
* Instead of using the CutFEM-characterization of the space, one could
  use an XFEM-characterization, cf. nxfem.py.
* To obtain higher order accuracy (also w.r.t. the geometry
  approximation) isoparametric unfitted methods can be used, cf.
  nxfem_higher_order.py.

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

from math import pi


# -------------------------------- PARAMETERS ---------------------------------
# Domain corners
ll, ur = (-1.5, -1.5), (1.5, 1.5)
# mesh size
maxh = 0.2
# Finite element order
order = 2

# Diffusion coefficients for the sub-domains (NEG/POS):
alpha = [1.0, 2.0]
# Nitsche penalty parameter
lambda_nitsche = 20


# ----------------------------------- MAIN ------------------------------------
# We generate the background mesh of the domain and use a simplicity
# triangulation to obtain a mesh with quadrilaterals use
# 'quad_dominated=True'
square = SplineGeometry()
square.AddRectangle(ll, ur, bc=1)
mesh = Mesh(square.GenerateMesh(maxh=maxh, quad_dominated=False))

# Manufactured solution and corresponding r.h.s. data CoefficientFunctions:
r44 = x**4 + y**4
r41 = sqrt(sqrt(x**4 + y**4))
r4m3 = 1.0 / (r41 * r41 * r41)
r66 = x**6 + y**6
r63 = sqrt(r66)
r22 = x**2 + y**2
r21 = sqrt(r22)
solution = [1 + pi / 2 - sqrt(2.0) * cos(pi / 4 * r44), pi / 2 * r41]
coef_f = [(-1 * sqrt(2) * pi * (pi * cos(pi / 4 * (r44)) * (r66)
                                + 3 * sin(pi / 4 * (r44)) * (r22))),
          (-2 * pi * 3 / 2 * (r4m3) * (-(r66) / (r44) + (r22)))]


# Level set function of the domain (phi = ||x||_4 - 1) and its
# interpolation:
levelset = sqrt(sqrt(x**4 + y**4)) - 1.0

if order > 1:
    from xfem.lsetcurv import LevelSetMeshAdaptation
    lsetadap = LevelSetMeshAdaptation(mesh, order=order, levelset=levelset)
else:
    lsetadap = NoDeformation(mesh, levelset)
lsetp1 = lsetadap.lset_p1

# Background FESpaces (used as CutFESpaces later-on):
Vh = H1(mesh, order=order, dirichlet=[1, 2, 3, 4])
VhG = FESpace([Vh, Vh])
print("unknowns in background FESpace (2 x standard unknowns): ", VhG.ndof)

# Gathering information on cut elements:
#  * domain of (volume/boundary) element:
#    * NEG: only negative level set values
#    * POS: only positive level set values
#    * IF: cut element (negative and positive) level set values
#  * cut ratio:
#    If element is cut this describes the ratio between the measure of
#    part in the negative domain and the measure of the full element.
ci = CutInfo(mesh, lsetp1)

# Overwrite freedofs (degrees of freedoms that should be solved for) of
# VhG to mark only dofs that are involved in the cut problem. Use cut
# information of ci here:
# (Note: "hasneg": has (also) negative level set values and
#        "haspos": has (also) positive level set values)
hasneg = ci.GetElementsOfType(HASNEG)
haspos = ci.GetElementsOfType(HASPOS)
freedofs = VhG.FreeDofs()
freedofs &= CompoundBitArray([GetDofsOfElements(Vh, hasneg),
                              GetDofsOfElements(Vh, haspos)])

# coefficients / parameters:
n = 1.0 / grad(lsetp1).Norm() * grad(lsetp1)
h = specialcf.mesh_size
# the cut ratio extracted from the cutinfo-class
kappa = (CutRatioGF(ci), 1.0 - CutRatioGF(ci))
# Nitsche stabilization parameter:
stab = lambda_nitsche * (alpha[1] + alpha[0]) / h

# expressions of test and trial functions (u and v are tuples):
u = VhG.TrialFunction()
v = VhG.TestFunction()

gradu = [grad(ui) for ui in u]
gradv = [grad(vi) for vi in v]

average_flux_u = sum([- kappa[i] * alpha[i] * gradu[i] * n for i in [0, 1]])
average_flux_v = sum([- kappa[i] * alpha[i] * gradv[i] * n for i in [0, 1]])

# Integration domains for integration on negative/positive sub-domains
# and on the interface: Here, the integration is (geometrically) exact
# if the "levelset"-argument is a piecewise (multi-)linear function. The
# integration order is chosen according to the arguments in the
# multi-linear forms (but can be overwritten with "force_intorder" in
# the integration domain). If the "levelset"-argument is not a
# (multi-)linear function, you can use the "subdivlvl" argument to add
# additional refinement levels for the geometry approximation.
lset_neg = {"levelset": lsetp1, "domain_type": NEG, "subdivlvl": 0}
lset_pos = {"levelset": lsetp1, "domain_type": POS, "subdivlvl": 0}
lset_if = {"levelset": lsetp1, "domain_type": IF, "subdivlvl": 0}

# bilinear forms:
a = BilinearForm(VhG, symmetric=True)
# l.h.s. domain integrals:
a += SymbolicBFI(levelset_domain=lset_neg, form=alpha[0] * gradu[0] * gradv[0])
a += SymbolicBFI(levelset_domain=lset_pos, form=alpha[1] * gradu[1] * gradv[1])
# Nitsche integrals:
a += SymbolicBFI(levelset_domain=lset_if,
                 form=(average_flux_u * (v[0] - v[1]) +
                       average_flux_v * (u[0] - u[1]) +
                       stab * (u[0] - u[1]) * (v[0] - v[1])))

f = LinearForm(VhG)
# r.h.s. domain integrals:
f += SymbolicLFI(levelset_domain=lset_neg, form=coef_f[0] * v[0])
f += SymbolicLFI(levelset_domain=lset_pos, form=coef_f[1] * v[1])

# solution vector
gfu = GridFunction(VhG)

with lsetadap:
    # setting domain boundary conditions:
    gfu.components[1].Set(solution[1], BND)

    # setting up matrix and vector
    a.Assemble()
    f.Assemble()

# homogenization of boundary data and solution of linear system
rhs = gfu.vec.CreateVector()
rhs.data = f.vec - a.mat * gfu.vec
update = gfu.vec.CreateVector()
update.data = a.mat.Inverse(freedofs) * rhs
gfu.vec.data += update

# visualization of (discrete) solution: Wherever (interpolated) level
# set function is negative visualize the first component, where it is#
# positive visualize the second component
u_coef = IfPos(lsetp1, gfu.components[1], gfu.components[0])

# visualize levelset, interpolated levelset and discrete solution:
# (Note that the visualization does not respect the discontinuities.
# They are smeared out. To see kinks or jumps more clearly increase the
# subdivision option in the visualization.)
Draw(levelset, mesh, "levelset")
Draw(lsetp1, mesh, "levelset_P1")
Draw(u_coef, mesh, "u")

# Error coefficients:
err_sqr_coefs = [(gfu.components[i] - solution[i])**2 for i in [0, 1]]

# Computation of L2 error:
with lsetadap:
    l2error = sqrt(Integrate(levelset_domain=lset_neg, cf=err_sqr_coefs[0],
                             mesh=mesh, order=2 * order)
                   + Integrate(levelset_domain=lset_pos, cf=err_sqr_coefs[1],
                               mesh=mesh, order=2 * order))

print("L2 error : ", l2error)
