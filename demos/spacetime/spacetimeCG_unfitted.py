"""
In this example we solve a scalar *unfitted* PDE problem on a moving
domain. The setting is the same as in `spacetimeDG_unfitted.py`
However, here, we apply a continuous-in-time Petrov-Galerkin method.

Domain + PDE problem:
---------------------
As in `spacetimeDG_unfitted.py`

Discretisation:
---------------
* Background space-time finite element space restricted to active domain
* Ghost penalty stabilization to deal with bad cuts (version as in [1])
  and in order to define a proper extension to a neighborhood, cf. [2]

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
from netgen.geom2d import SplineGeometry
from xfem import *
from math import pi
from xfem.lset_spacetime import *
ngsglobals.msg_level = 1

# -------------------------------- PARAMETERS ---------------------------------

# DISCRETIZATION PARAMETERS:
# parameter for refinement study:
i = 2
n_steps = 2**i
space_refs = i

# polynomial order in time
k_t = 3
# polynomial order in space
k_s = k_t
# polynomial order in time for level set approximation
lset_order_time = k_t
# integration order in time
time_order = 2 * k_t
# time stepping parameters
tstart = 0
tend = 0.5
delta_t = (tend - tstart) / n_steps
maxh = 0.5
# ghost penalty parameter
gamma = 0.05
# map from reference time to physical time
told = Parameter(tstart)
t = told + delta_t * tref

# PROBLEM SETUP:

# outer domain:
rect = SplineGeometry()
rect.AddRectangle([-0.6, -1], [0.6, 1])

# level set geometry
# radius of disk (the geometry)
R = 0.5
# position shift of the geometry in time
rho = (1 / (pi)) * sin(2 * pi * t)
# convection velocity:
w = CoefficientFunction((0, rho.Diff(t)))
max_velocity = 2
# level set
r = sqrt(x**2 + (y - rho)**2)
levelset = r - R

# diffusion coeff
alpha = 1
# solution
u_exact = cos(pi * r / R) * sin(pi * t)
# r.h.s.
coeff_f = (u_exact.Diff(t)
           - alpha * (u_exact.Diff(x).Diff(x) + u_exact.Diff(y).Diff(y))
           + w[0] * u_exact.Diff(x) + w[1] * u_exact.Diff(y)).Compile()

# ----------------------------------- MAIN ------------------------------------

ngmesh = rect.GenerateMesh(maxh=maxh, quad_dominated=False)
for j in range(space_refs):
    ngmesh.Refine()
mesh = Mesh(ngmesh)

# spatial FESpace for solution
fes = H1(mesh, order=k_s, dgjumps=True)
# time finite elements (nodal!)
tfe_i = ScalarTimeFE(k_t, skip_first_node=True)  # interior shapes
tfe_e = ScalarTimeFE(k_t, only_first_node=True)  # exterior shapes (init. val.)
tfe_t = ScalarTimeFE(k_t - 1)                    # test shapes
# space-time finite element space
st_fes_i, st_fes_e, st_fes_t = [tfe * fes for tfe in [tfe_i, tfe_e, tfe_t]]

# Space time version of Levelset Mesh Adaptation object. Also offers integrator
# helper functions that involve the correct mesh deformation
lsetadap = LevelSetMeshAdaptation_Spacetime(mesh, order_space=k_s,
                                            order_time=lset_order_time,
                                            threshold=0.5,
                                            discontinuous_qn=True)

# lset epsilon perturbation for extended facet patch bfi domain
eps = 2 * max_velocity * delta_t

gfu_i = GridFunction(st_fes_i)
gfu_e = GridFunction(st_fes_e)

u_last = CreateTimeRestrictedGF(gfu_e, 0)

u_i = st_fes_i.TrialFunction()
u_e = st_fes_e.TrialFunction()
v_t = st_fes_t.TestFunction()

scene = DrawDC(lsetadap.levelsetp1[TOP], u_last, 0, mesh, "u_last",
               deformation=lsetadap.deformation[TOP])

lset_p1_slice = GridFunction(lsetadap.levelsetp1[BOTTOM].space)

h = specialcf.mesh_size

ba_strip = BitArray(mesh.ne)
ba_facets = BitArray(mesh.nfacet)
ci = CutInfo(mesh, time_order=time_order)
ci_slice = CutInfo(mesh)

dQ = delta_t * dCut(lsetadap.levelsetp1[INTERVAL], NEG, time_order=2 * k_t,
                    deformation=lsetadap.deformation[INTERVAL],
                    definedonelements=ci.GetElementsOfType(HASNEG))
dOmnew = dCut(lsetadap.levelsetp1[TOP], NEG,
              deformation=lsetadap.deformation[TOP],
              definedonelements=ci.GetElementsOfType(HASNEG))
dw = delta_t * dFacetPatch(definedonelements=ba_facets, time_order=time_order,
                           deformation=lsetadap.deformation[INTERVAL])


def dt(u): return 1.0 / delta_t * dtref(u)


a_i = BilinearForm(trialspace=st_fes_i, testspace=st_fes_t, check_unused=False)

a_i += v_t * (dt(u_i) - dt(lsetadap.deform) * grad(u_i)) * dQ
a_i += (alpha * InnerProduct(grad(u_i), grad(v_t))) * dQ
a_i += (v_t * InnerProduct(w, grad(u_i))) * dQ
a_i += h**(-2) * (1 + delta_t / h) * gamma * \
    (u_i - u_i.Other()) * (v_t - v_t.Other()) * dw

f = LinearForm(st_fes_t)
f += coeff_f * v_t * dQ
f += -v_t * (dt(gfu_e) - dt(lsetadap.deform) * grad(gfu_e)) * dQ
f += -(alpha * InnerProduct(grad(gfu_e), grad(v_t))) * dQ
f += -(v_t * InnerProduct(w, grad(gfu_e))) * dQ

# set initial values
u_last.Set(fix_tref(u_exact, 0))
# project u_last at the beginning of each time step
lsetadap.ProjectOnUpdate(u_last)

ba_plus_hasneg_old, els_test = BitArray(mesh.ne), BitArray(mesh.ne)
ba_plus_hasneg_old.Set()

while tend - told.Get() > delta_t / 2:
    lsetadap.CalcDeformation(levelset)
    gfu_e.Set(u_last)

    # update markers in (space-time) mesh
    ci.Update(lsetadap.levelsetp1[INTERVAL], time_order=time_order)

    # re-evaluate the "active dofs" in the space time slab
    InterpolateToP1(lsetadap.levelsetp1[BOTTOM] - eps, lset_p1_slice)
    ci_slice.Update(lset_p1_slice)
    ba_plus_hasneg = BitArray(ci_slice.GetElementsOfType(HASNEG))

    InterpolateToP1(lsetadap.levelsetp1[BOTTOM] + eps, lset_p1_slice)
    ci_slice.Update(lset_p1_slice)
    ba_minus_haspos = BitArray(ci_slice.GetElementsOfType(HASPOS))

    ba_strip[:] = ba_minus_haspos & ba_plus_hasneg
    ba_facets[:] = GetFacetsWithNeighborTypes(mesh, a=ba_strip,
                                              b=ba_plus_hasneg)
    active_dofs = GetDofsOfElements(st_fes_i, ba_plus_hasneg)

    # Check element history for method of lines time-derivative approx.
    els_test[:] = ci.GetElementsOfType(HASNEG) & ~ba_plus_hasneg_old
    assert sum(els_test) == 0, 'Some active elements do not have a history.\
        You might want to increase eps'

    ba_plus_hasneg_old[:] = ba_plus_hasneg

    a_i.Assemble()
    f.Assemble()

    # solve linear system
    gfu_i.vec.data = a_i.mat.Inverse(active_dofs) * f.vec

    # evaluate upper trace of solution for
    #  * for error evaluation
    #  * upwind-coupling to next time slab
    RestrictGFInTime(spacetime_gf=gfu_i, reference_time=1.0, space_gf=u_last)

    # compute error at final time
    l2error = sqrt(
        Integrate((fix_tref(u_exact, 1) - u_last)**2 * dOmnew, mesh))

    # update time variable (ParameterCL)
    told.Set(told.Get() + delta_t)
    print("\rt = {0:12.9f}, L2 error = {1:12.9e}".format(told.Get(), l2error))

    try:
        __builtin__
        __IPYTHON__
        scene.Redraw()
    except NameError:
        scene.Redraw(blocking=True)
