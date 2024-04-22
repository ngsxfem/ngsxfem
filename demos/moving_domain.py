"""
In this example we consider a scalar time-dependent convection-diffusion
problem on a moving domain. The discretisation uses a level set based
geometry description, a cut finite element method for the spatial
discretisation and utilised ghost-penalty stabilisation to realise a
discrete extension of the solution into an extension strip. With this
extension, a method of lines approximation of the time-derivative is
possible to create an Eulerian time-stepping method. For the numerical
analysis of this method, see [1]. In contrast to [1], this file uses
an isoparametric mapping approach to realise a higher-order geometry
approximation.

Used Features:
--------------
* Level set geometry description and element markings, see the 'basic'
  jupyter tutorial.

* Integrators with respect to level set geometry, see the 'cutfem'
  jupyter tutorial.

* Higher-order geometry approximation, see the 'cutfem' jupyter
  tutorial.

* Ghost-penalty stabilisation, see the 'spacetime' jupyter tutorial.

Literature:
-----------
[1] C. Lehrenfeld, M. Olshanskii, An Eulerian finite element method for
    PDEs in time-dependent domains. ESAIM M2AN 53(2):585-614, 2019.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.lsetcurv import *

from math import pi, ceil

ngsglobals.msg_level = 2
SetNumThreads(4)

# -------------------------------- PARAMETERS ---------------------------------
# Finite element space order
k = 1
# Initial mesh size and time-step size
h0, t0 = 0.2, 0.1
# Number of mesh and time-step bisections
Lx, Lt = 3, 3

# Background geometry corners
lowerleft, upperright = [-0.7, -0.7], [0.9, 0.7]

# Time parameters
dt = t0 * 0.5**Lt
T_end = 0.4
t = Parameter(0.0)

# Problem and discretisation parameters
nu = 1e-5
c_gamma = 1

# Sparse direct solver
inverse = ""

# Problem Data
# Initial condition and right-hand side
rho = CoefficientFunction(1 / pi * sin(2 * pi * t))

r0 = 0.5
r1 = pi / (2 * r0)

u0 = CoefficientFunction(cos(r1 * sqrt(x**2 + y**2))**2)
u_ex = CoefficientFunction(cos(r1 * sqrt((x - rho)**2 + y**2))**2)

grad_u_ex = CoefficientFunction((-pi * sin(pi / r0 * sqrt((x - rho)**2 + y**2))
                                 * (x - rho) / sqrt((x - rho)**2 + y**2),
                                 - pi * sin(pi / r0 *
                                            sqrt((x - rho)**2 + y**2))
                                 * y / sqrt((x - rho)**2 + y**2)))
# rhs from py_tutorial
rhs = nu * CoefficientFunction(-(pi / r0) * r1
                               * (sin(r1 * sqrt((x - rho)**2 + y**2))
                                  * sin(r1 * sqrt((x - rho)**2 + y**2))
                                  - cos(r1 * sqrt((x - rho)**2 + y**2))
                                  * cos(r1 * sqrt((x - rho)**2 + y**2)))
                               + (pi / r0)
                               * cos(r1 * sqrt((x - rho)**2 + y**2))
                               * sin(r1 * sqrt((x - rho)**2 + y**2))
                               * (1 / sqrt((x - rho)**2 + y**2)))

# Convection Field
w = CoefficientFunction((2 * cos(2 * pi * t), 0))
div_w = 0.0
velmax = 2


def Rho(t):
    return CoefficientFunction(1 / pi * sin(2 * pi * t))


# Level set
def levelset_func(t):
    return sqrt((x - Rho(t))**2 + y**2) - r0


# ----------------------------------- MAIN ------------------------------------
background_domain = SplineGeometry()
background_domain.AddRectangle(lowerleft, upperright, bc=1)
ngmesh = background_domain.GenerateMesh(maxh=h0, quad_dominated=False)
for i in range(Lx):
    ngmesh.Refine()
mesh = Mesh(ngmesh)
h_max = h0 * 0.5**Lx

# FE Space
V = H1(mesh, order=k, dgjumps=True)
gfu = GridFunction(V)


# Higher order discrete level set approximation
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=k, threshold=0.1,
                                      discontinuous_qn=True)
deformation = lsetmeshadap.deform
lsetp1 = lsetmeshadap.lset_p1

# Cut-Info classes for element marking
ci_main = CutInfo(mesh)
ci_inner = CutInfo(mesh)
ci_outer = CutInfo(mesh)

# Element and facet markers
els_hasneg = ci_main.GetElementsOfType(HASNEG)
els_outer = ci_outer.GetElementsOfType(HASNEG)
els_inner = ci_inner.GetElementsOfType(NEG)

els_ring = BitArray(mesh.ne)
facets_ring = BitArray(mesh.nedge)
els_outer_old, els_test = BitArray(mesh.ne), BitArray(mesh.ne)

# Discretisation
delta = 1 * dt * velmax
K_tilde = int(ceil(delta / h_max))
gamma_s = c_gamma * K_tilde

u, v = V.TnT()
h = specialcf.mesh_size

dx = dCut(levelset=lsetp1, domain_type=NEG, definedonelements=els_hasneg,
          deformation=deformation)
dw = dFacetPatch(definedonelements=facets_ring, deformation=deformation)

# Bilinear and linear forms
a = RestrictedBilinearForm(V, element_restriction=els_outer,
                           facet_restriction=facets_ring, check_unused=False)
a += (1 / dt) * u * v * dx
a += nu * InnerProduct(Grad(u), Grad(v)) * dx
a += (InnerProduct(w, grad(u)) * v + div_w * u * v) * dx
a += gamma_s * (1 / h**2) * (u - u.Other()) * (v - v.Other()) * dw

f = LinearForm(V)
f += rhs * v * dx
f += (1 / dt) * gfu * v * dx


# Error computation
errors_L2, errors_H1 = [], []
dx_2k = dx.order(2 * k)


def CompErrs():
    l2 = sqrt(Integrate((gfu - u_ex)**2 * dx_2k, mesh))
    h1 = sqrt(Integrate((grad(gfu) - grad_u_ex)**2 * dx_2k, mesh))

    errors_L2.append(l2)
    errors_H1.append(h1)

    return l2


# Project the solution defined with respect to the last mesh deformation
# onto the the mesh with the current mesh deformation.
lsetmeshadap.ProjectOnUpdate(gfu)

# Time stepping loop
gfu.Set(u0)
els_outer_old.Set()

Draw(gfu, mesh, "gfu")
Draw(gfu - u_ex, mesh, "l2err")

for it in range(1, int(T_end / dt + 0.5) + 1):
    t.Set(it * dt)
    # Update mesh deformation (Note: this updates gfu)
    deformation = lsetmeshadap.CalcDeformation(levelset_func(t))

    # Update Element markers: extension facets. Updating the cut-info
    # classes automatically updates the marker BitArrays.
    InterpolateToP1(levelset_func(t) - delta, lsetp1)
    ci_outer.Update(lsetp1)
    InterpolateToP1(levelset_func(t) + delta, lsetp1)
    ci_inner.Update(lsetp1)

    els_ring.Clear(), facets_ring.Clear()
    els_ring |= els_outer & ~els_inner
    facets_ring |= GetFacetsWithNeighborTypes(mesh, a=els_outer, b=els_ring)

    # Update Element markers: Domain elements
    InterpolateToP1(levelset_func(t), lsetp1)
    ci_main.Update(lsetp1)

    active_dofs = GetDofsOfElements(V, els_outer)

    # Check element history for method of lines time-derivative approx.
    els_test[:] = els_hasneg & ~els_outer_old
    assert sum(els_test) == 0, 'Some active elements do not have a history'

    els_outer_old[:] = els_outer

    # Update Linear System
    a.Assemble(reallocate=True)
    f.Assemble()
    inv = a.mat.Inverse(active_dofs, inverse=inverse)
    # Solve Problem
    gfu.vec.data = inv * f.vec

    # Compute Errors
    l2err = CompErrs()

    print(f"Lx = {Lx}, dt = {dt:8.6f}, t = {t.Get():8.6f}, ", end="")
    print(f"error(L2) = {l2err:4.2e}, ", end="")
    print(f"active_els = {sum(els_outer)}, K = {K_tilde}")
    Redraw(blocking=True)


# Post Processing
err_l2l2 = sqrt(dt * sum([e**2 for e in errors_L2]))
err_l2h1 = sqrt(dt * sum([e**2 for e in errors_H1]))
# err_l1l2 = sum([e*dt for e in errors_L2])
# err_l1h1 = sum([e*dt for e in errors_H1])

print("\n--------------------------------------------------------")
print("  L2(0,T;L2) Error = {:6.4e}".format(err_l2l2))
print("Linf(0,T;L2) Error = {:6.4e}".format(max(errors_L2)))
print("  L2(0,T;H1) Error = {:6.4e}".format(err_l2h1))
