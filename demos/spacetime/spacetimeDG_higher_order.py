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

## DISCRETIZATION PARAMETERS:
# parameter for refinement study:
i=2
n_steps = 2**i
space_refs = i

# polynomial order in time
k_t = 2
# polynomial order in space
k_s = k_t
# polynomial order in time for level set approximation
lset_order_time = k_t
# integration order in time
time_order = 2*k_t
# time stepping parameters
tstart = 0
tend = 0.5
delta_t = (tend - tstart)/n_steps
maxh = 0.5
# ghost penalty parameter
gamma = 0.05
# map from reference time to physical time
told = Parameter(tstart)
t = told + delta_t*tref

## PROBLEM SETUP:

# outer domain:
rect = SplineGeometry()
rect.AddRectangle([-0.6,-1],[0.6,1])

# level set geometry 
# radius of disk (the geometry)
R = 0.5
# position shift of the geometry in time
rho =  (1/(pi))*sin(2*pi*t)
#convection velocity:
w = CoefficientFunction((0,rho.Diff(t)))
# level set
r = sqrt(x**2+(y-rho)**2)
levelset= r - R

# diffusion coeff
alpha = 1
# solution
u_exact = cos(pi*r/R) * sin(pi*t)
# r.h.s.
coeff_f = (u_exact.Diff(t) - alpha * (u_exact.Diff(x).Diff(x)+u_exact.Diff(y).Diff(y)) \
           + w[0]*u_exact.Diff(x) + w[1]*u_exact.Diff(y)).Compile()

# ----------------------------------- MAIN ------------------------------------

ngmesh = rect.GenerateMesh(maxh=maxh, quad_dominated=False)
for j in range(space_refs):
    ngmesh.Refine()
mesh = Mesh (ngmesh)

# spatial FESpace for solution
fes1 = H1(mesh, order=k_s)
# time finite element (nodal!)
tfe = ScalarTimeFE(k_t) 
# space-time finite element space
st_fes = SpaceTimeFESpace(fes1,tfe, flags = {"dgjumps": True})

# Space time version of Levelset Mesh Adapation object. Also offers integrator helper functions that
# involve the correct mesh deformation 
lsetadap = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = lset_order_time,
                                            threshold=0.5, discontinuous_qn=True)

gfu = GridFunction(st_fes)
u_last = CreateTimeRestrictedGF(gfu,0)

Draw(IfPos(lsetadap.levelsetp1[TOP],0,u_last),mesh,"u_last")

u,v = st_fes.TnT()
h = specialcf.mesh_size

hasneg_ints_a = [lsetadap.BFI(NEG,INTERVAL, v*(dt(u) + lsetadap.mesh_velocity*grad(u) )),
                 lsetadap.BFI(NEG,INTERVAL, alpha * delta_t* InnerProduct( grad(u), grad(v)) ),
                 lsetadap.BFI(NEG,INTERVAL, delta_t*v*InnerProduct(w,grad(u))),
                 lsetadap.BFI(NEG,BOTTOM, fix_t(u,0)*fix_t(v,0))]
patch_ints_a = [SymbolicFacetPatchBFI(form = h**(-2)*delta_t*(1+delta_t/h)*gamma*(u - u.Other() )*(v - v.Other()),
                                      skeleton=False, time_order=time_order, deformation=lsetadap.deform)]
hasneg_ints_f = [lsetadap.LFI(NEG,INTERVAL, delta_t*coeff_f*v),
                 lsetadap.LFI(NEG,BOTTOM, u_last*fix_t(v,0))]

a = BilinearForm(st_fes,"a", check_unused=False)
for integrator in hasneg_ints_a + patch_ints_a:
    a += integrator
f = LinearForm(st_fes)
for integrator in hasneg_ints_f:
    f += integrator

# set initial values
u_last.Set(fix_t(u_exact,0))
# project u_last at the beginning of each time step
lsetadap.ProjectOnUpdate(u_last)
ci = CutInfo(mesh,time_order=time_order)

while tend - told.Get() > delta_t/2:
    lsetadap.CalcDeformation(levelset)

    # update markers in (space-time) mesh
    ci.Update(lsetadap.levelsetp1[INTERVAL],time_order=time_order)

    # re-compute the facets for stabilization:
    ba_facets = GetFacetsWithNeighborTypes(mesh,a=ci.GetElementsOfType(HASNEG),
                                                b=ci.GetElementsOfType(IF))
    active_dofs = GetDofsOfElements(st_fes,ci.GetElementsOfType(HASNEG))

    # re-set definedonelements-markers according to new markings:
    for integrator in hasneg_ints_a + hasneg_ints_f:
        integrator.SetDefinedOnElements(ci.GetElementsOfType(HASNEG))
    for integrator in patch_ints_a:
        integrator.SetDefinedOnElements(ba_facets)

    a.Assemble()
    f.Assemble()
    
    # solve linear system
    inv = a.mat.Inverse(active_dofs)
    gfu.vec.data = inv * f.vec.data
    
    # evaluate upper trace of solution for
    #  * for error evaluation 
    #  * upwind-coupling to next time slab
    RestrictGFInTime(spacetime_gf=gfu,reference_time=1.0,space_gf=u_last)
    
    # compute error at final time
    l2error = sqrt(lsetadap.Integrate(NEG,TOP,(fix_t(u_exact,1) - u_last)**2, order=2*k_s))

    # update time variable (ParameterCL)
    told.Set(told.Get() + delta_t)
    print("\rt = {0:12.9f}, L2 error = {1:12.9e}".format(told.Get(),l2error))
    
    Redraw(blocking=True)

