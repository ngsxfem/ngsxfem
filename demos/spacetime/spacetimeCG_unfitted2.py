"""
unfitted Heat equation with Neumann b.c. solved with an unfitted isoparametric
space-time discretisation.
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
i = 3
n_steps = 2**i
space_refs = 0

# polynomial order in time
k_t = 2
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
fes1 = H1(mesh, order=k_s)
# time finite elements (nodal!)
tfe = ScalarTimeFE(k_t) 
tfe_i = ScalarTimeFE(k_t, skip_first_node=True) # interior
tfe_e = ScalarTimeFE(k_t, only_first_node=True) # exterior (inital values)
tfe_t = ScalarTimeFE(k_t-1)                     # test
# space-time finite element space
st_fes = SpaceTimeFESpace(fes1,tfe, flags = {"dgjumps": True})
st_fes_i = SpaceTimeFESpace(fes1,tfe_i, flags = {"dgjumps": True})
st_fes_e = SpaceTimeFESpace(fes1,tfe_e, flags = {"dgjumps": True})
st_fes_t = SpaceTimeFESpace(fes1,tfe_t, flags = {"dgjumps": True})

# Space time version of Levelset Mesh Adapation object. Also offers integrator
# helper functions that involve the correct mesh deformation
lsetadap = LevelSetMeshAdaptation_Spacetime(mesh, order_space=k_s,
                                            order_time=lset_order_time,
                                            threshold=0.5,
                                            discontinuous_qn=True)

## lset epsilon pertubation for extended facet patch bfi domain
eps=0.1*delta_t

gfu_i = GridFunction(st_fes_i)
gfu_e = GridFunction(st_fes_e)

u_last = CreateTimeRestrictedGF(gfu_e,0)
#u_last_shifted = GridFunction(fes1)

u_i = st_fes_i.TrialFunction()
u_e = st_fes_e.TrialFunction()
v_t = st_fes_t.TestFunction()

scene = DrawDC(lsetadap.levelsetp1[TOP], u_last, 0, mesh, "u_last",
               deformation=lsetadap.deformation[TOP])

lset_p1_slice = GridFunction(fes1)

h = specialcf.mesh_size

ba_facets = BitArray(mesh.nfacet)
ci = CutInfo(mesh, time_order=time_order)
ci_slice = CutInfo(mesh)

dQ = delta_t * dCut(lsetadap.levelsetp1[INTERVAL], NEG, time_order=2 * k_t,
                    deformation=lsetadap.deformation[INTERVAL],
                    definedonelements=ci.GetElementsOfType(HASNEG))
#dOmold = dCut(lsetadap.levelsetp1[BOTTOM], NEG,
              #deformation=lsetadap.deformation[BOTTOM],
              #definedonelements=ci.GetElementsOfType(HASNEG))
dOmnew = dCut(lsetadap.levelsetp1[TOP], NEG,
              deformation=lsetadap.deformation[TOP],
              definedonelements=ci.GetElementsOfType(HASNEG))
dw = delta_t * dFacetPatch(definedonelements=ba_facets, time_order=time_order,
                           deformation=lsetadap.deformation[INTERVAL])


def dt(u): return 1.0 / delta_t * dtref(u)

a_i = BilinearForm(trialspace = st_fes_i, testspace = st_fes_t, check_unused=False, symmetric=False)
a_e = BilinearForm(trialspace = st_fes_e, testspace = st_fes_t, check_unused=False, symmetric=False)

a_i += v_t * (dt(u_i) - dt(lsetadap.deform) * grad(u_i)) * dQ
a_i += (alpha * InnerProduct(grad(u_i), grad(v_t))) * dQ
a_i += (v_t * InnerProduct(w, grad(u_i))) * dQ
#a += (fix_tref(u, 0) * fix_tref(v, 0)) * dOmold
a_i += h**(-2) * (1 + delta_t / h) * gamma * \
    (u_i - u_i.Other()) * (v_t - v_t.Other()) * dw

a_e += v_t * (dt(u_e) - dt(lsetadap.deform) * grad(u_e)) * dQ
a_e += (alpha * InnerProduct(grad(u_e), grad(v_t))) * dQ
a_e += (v_t * InnerProduct(w, grad(u_e))) * dQ

f = LinearForm(st_fes_t)
f += coeff_f * v_t * dQ
#f += u_last * fix_tref(v, 0) * dOmold

# set initial values
u_last.Set(fix_tref(u_exact, 0))
# project u_last at the beginning of each time step
lsetadap.ProjectOnUpdate(u_last)

while tend - told.Get() > delta_t / 2:
    lsetadap.CalcDeformation(levelset)

    # update markers in (space-time) mesh
    ci.Update(lsetadap.levelsetp1[INTERVAL], time_order=time_order)
    
    # re-evaluate the "active dofs" in the space time slab
    active_dofs = GetDofsOfElements(st_fes,ci.GetElementsOfType(HASNEG))
        
    # Actually, one could also subtract epsilon from all vector elements...
    InterpolateToP1(lsetadap.levelsetp1[BOTTOM]-eps,lset_p1_slice)
    ci_slice.Update(lset_p1_slice)
    ba_plus_hasneg = BitArray(ci_slice.GetElementsOfType(HASNEG))
        
    InterpolateToP1(lsetadap.levelsetp1[BOTTOM]+eps,lset_p1_slice)
    ci_slice.Update(lset_p1_slice)
    ba_minus_haspos = BitArray(ci_slice.GetElementsOfType(HASPOS))
        
    ba_strip = BitArray(ba_minus_haspos & ba_plus_hasneg)
    ba_facets[:] = GetFacetsWithNeighborTypes(mesh,a=ba_strip,b=ba_plus_hasneg)
    active_dofs |= GetDofsOfElements(st_fes,ba_strip)

    #print(active_dofs)

    a_i.Assemble(reallocate=True)
    a_e.Assemble(reallocate=True)
    f.Assemble()

    #print("a_e: ", a_e.mat)

    # solve linear system
    inv = a_i.mat.Inverse(active_dofs)
    f.vec.data -= a_e.mat * gfu_e.vec
    gfu_i.vec.data =  inv * f.vec
    
    # evaluate upper trace of solution for
    #  * for error evaluation
    #  * upwind-coupling to next time slab
    RestrictGFInTime(spacetime_gf=gfu_i, reference_time=1.0, space_gf=u_last)
    gfu_e.Set(u_last)

    # compute error at final time
    l2error = sqrt(
        Integrate((fix_tref(u_exact, 1) - u_last)**2 * dOmnew, mesh))

    # update time variable (ParameterCL)
    told.Set(told.Get() + delta_t)
    print("\rt = {0:12.9f}, L2 error = {1:12.9e}".format(told.Get(), l2error))
    
    #break
    
    try:
        __builtin__
        __IPYTHON__
        scene.Redraw()
    except NameError:
        scene.Redraw(blocking=True)
