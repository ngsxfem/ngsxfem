"""
In this example we solve a scalar *unfitted* PDE problem on a moving
domain. A discontinuous-in-time space-time formulation is applied.
Natural boundary conditions are applied which simplifies the variational
formulation. To stabilize arbitrary cut configurations, we use a
space-time version of the ghost penalty method.

Domain:
-------
The background domain is [-0.6,0.6]x[-1,1]x[0,0.5] (2D+time interval)
while the physical domain is a circle that is traveling up and down
over time.

PDE problem:
------------
  u_t + wx·u_x + wy·u_y - (u_xx + u_yy) = f in  Omega(t) (where lset is neg.)
        (-u_x+wx u)·nx+(-u_y+wy u)·ny u = 0 on dOmega(t) (where lset is zero.)
where w = (wx,wy) is a divergence-free vector field.
The r.h.s. term f is chosen according to a manufactured solution.

Discretisation:
---------------
* Background space-time finite element space restricted to active domain

* Ghost penalty stabilization to deal with bad cuts (version as in [1])

Implementational aspects (cf. [1] and [2] for details):
-------------------------------------------------------
* Geometry approximation in space-time using isoparametric unfitted FEM

* Projection operator that maps solutions from one deformed mesh to another

* A (sparse) direct solver is applied to solve the arising linear systems.

References:
-----------
All concepts that are used here are explained in the jupyter-tuorials
`spacetime.ipynb`. As a simplified setting without cut configurations,
we also refer to the `spacetimeDG_fitted.py` demo.

Literature:
-----------
[1] J. Preuß, Higher order unfitted isoparametric space-time FEM on moving
    domains. Master's thesis, NAM, University of Göttingen, 2018.
[2] F. Heimann. On Discontinuous- and Continuous-In-Time Unfitted Space-Time
    Methods for PDEs on Moving Domains. Master's thesis, NAM, University of
    Göttingen, 2020.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from ngsolve import *
from ngsolve.meshes import OrthoBrick, Pnt, CSGeometry
from xfem import *
from math import pi
from xfem.lset_spacetime import *
ngsglobals.msg_level = 1

# -------------------------------- PARAMETERS ---------------------------------
# DISCRETIZATION PARAMETERS:

# Parameter for refinement study:
i = 2
n_steps = 2**i
space_refs = i

# Polynomial order in time
k_t = 2
# Polynomial order in space
k_s = k_t
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
geometry = CSGeometry()
geometry.Add(OrthoBrick(Pnt(-0.6, -1, -0.6), Pnt(0.6, 1, 0.6)))

D = 3

(xmin, xmax) = (-0.6, 0.6)
(ymin, ymax) = (-0.85, 0.85)
(zmin, zmax) = (-0.6, 0.6)

# Level set geometry
# Radius of sphere (the geometry)
R = 0.5
# Position shift of the geometry in time
rho = (1 / (pi)) * sin(2 * pi * t)
# Convection velocity:
w = CoefficientFunction((0, rho.Diff(t), 0))
max_velocity = 2

# Level set
r = sqrt(x**2 + (y - rho)**2 + z**2)
levelset = r - R

# Diffusion coefficient
alpha = 1
# Solution
u_exact = cos(pi * r / R) * sin(pi * t)
# R.h.s.
coeff_f = (u_exact.Diff(t)
           - alpha * (u_exact.Diff(x).Diff(x) + u_exact.Diff(y).Diff(y)
                      + u_exact.Diff(z).Diff(z))
           + w[0] * u_exact.Diff(x) + w[1] * u_exact.Diff(y)
           + w[2] * u_exact.Diff(z)).Compile()

# ----------------------------------- MAIN ------------------------------------
ngmesh = geometry.GenerateMesh(maxh=maxh, quad_dominated=False)
for j in range(space_refs):
    ngmesh.Refine()
mesh = Mesh(ngmesh)

# Spatial FESpace for solution
fes1 = L2(mesh, order=k_s, dgjumps=True)
# Time finite element (nodal!)
tfe = ScalarTimeFE(k_t)
# (Tensor product) space-time finite element space
st_fes = tfe * fes1

# Space time version of Levelset Mesh Adapation object. Also offers integrator
# helper functions that involve the correct mesh deformation
lsetadap = LevelSetMeshAdaptation_Spacetime(mesh, order_space=k_s,
                                            order_time=lset_order_time,
                                            threshold=0.5,
                                            discontinuous_qn=True)

gfu = GridFunction(st_fes)
u_last = CreateTimeRestrictedGF(gfu, 1)

scene = DrawDC(lsetadap.levelsetp1[TOP], u_last, 0, mesh, "u_last",
               deformation=lsetadap.deformation[TOP])

u, v = st_fes.TnT()
h = specialcf.mesh_size

ba_facets = BitArray(mesh.nfacet)
ba_facets_inner = BitArray(mesh.nfacet)
ci = CutInfo(mesh, time_order=0)

dQ = delta_t * dCut(lsetadap.levelsetp1[INTERVAL], NEG, time_order=time_order,
                    deformation=lsetadap.deformation[INTERVAL],
                    definedonelements=ci.GetElementsOfType(HASNEG))
dQ_f = delta_t * dCut(lsetadap.levelsetp1[INTERVAL], NEG,
                      time_order=time_order,
                      deformation=lsetadap.deformation[INTERVAL],
                      definedonelements=ba_facets_inner, skeleton=True)
dOmold = dCut(lsetadap.levelsetp1[BOTTOM], NEG,
              deformation=lsetadap.deformation[BOTTOM],
              definedonelements=ci.GetElementsOfType(HASNEG), tref=0)
dOmnew = dCut(lsetadap.levelsetp1[TOP], NEG,
              deformation=lsetadap.deformation[TOP],
              definedonelements=ci.GetElementsOfType(HASNEG), tref=1)
dw = delta_t * dFacetPatch(definedonelements=ba_facets, time_order=time_order,
                           deformation=lsetadap.deformation[INTERVAL])


def dt(u):
    return 1.0 / delta_t * dtref(u)


jump_u = u-u.Other()
jump_v = v-v.Other()
n = specialcf.normal(D)
mean_dudn = 0.5*n * (grad(u)+grad(u).Other())
mean_dvdn = 0.5*n * (grad(v)+grad(v).Other())

a = RestrictedBilinearForm(st_fes, "a", check_unused=False,
                           element_restriction=ci.GetElementsOfType(HASNEG),
                           facet_restriction=ba_facets_inner)
a += v * (dt(u) - dt(lsetadap.deform) * grad(u)) * dQ
a += (alpha * InnerProduct(grad(u), grad(v))) * dQ
a += (v * InnerProduct(w, grad(u))) * dQ
a += u * v * dOmold

a += (-mean_dvdn*jump_u - mean_dudn*jump_v + 50*k_s**2/h*jump_u*jump_v) * dQ_f

a += h**(-2) * (1 + delta_t / h) * gamma * \
    (u - u.Other()) * (v - v.Other()) * dw

f = LinearForm(st_fes)
f += coeff_f * v * dQ
f += u_last * v * dOmold

# Set initial values
u_last.Set(fix_tref(u_exact, 0))
# Project u_last at the beginning of each time step
lsetadap.ProjectOnUpdate(u_last)

while tend - told.Get() > delta_t / 2:
    lsetadap.CalcDeformation(levelset)

    # Update markers in (space-time) mesh
    ci.Update(lsetadap.levelsetp1[INTERVAL], time_order=0)

    # re-compute the facets for stabilization:
    ba_facets[:] = GetFacetsWithNeighborTypes(mesh,
                                              a=ci.GetElementsOfType(HASNEG),
                                              b=ci.GetElementsOfType(IF))
    ba_facets_inner[:] = GetFacetsWithNeighborTypes(
                                        mesh,
                                        a=ci.GetElementsOfType(HASNEG),
                                        b=ci.GetElementsOfType(HASNEG))
    active_dofs = GetDofsOfElements(st_fes, ci.GetElementsOfType(HASNEG))

    a.Assemble(reallocate=True)
    f.Assemble()

    # Solve linear system
    inv = a.mat.Inverse(active_dofs, inverse="")
    gfu.vec.data = inv * f.vec.data

    # Evaluate upper trace of solution for
    #  * for error evaluation
    #  * upwind-coupling to next time slab
    RestrictGFInTime(spacetime_gf=gfu, reference_time=1.0, space_gf=u_last)

    # Compute error at final time
    l2error = sqrt(Integrate((u_exact - u_last)**2 * dOmnew, mesh))

    # Update time variable (ParameterCL)
    told.Set(told.Get() + delta_t)
    print("\rt = {0:12.9f}, L2 error = {1:12.9e}".format(told.Get(), l2error))

    try:
        __builtin__
        __IPYTHON__
        scene.Redraw()
    except NameError:
        scene.Redraw(blocking=True)
