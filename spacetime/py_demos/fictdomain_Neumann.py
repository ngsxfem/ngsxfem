# unfitted Heat equation with Neumann b.c.
from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from math import pi

from xfem.lset_spacetime import *

ngsglobals.msg_level = 1

square = SplineGeometry()
square.AddRectangle([-1,-0.75],[1,1.5],bc=1)
ngmesh = square.GenerateMesh(maxh=0.05, quad_dominated=False)
mesh = Mesh (ngmesh)

n_ref = 0
for i in range(n_ref):
    mesh.Refine()

# polynomial order in time
k_t = 1
# polynomial order in space
k_s = 2
# spatial FESpace for solution
fes1 = H1(mesh, order=k_s, dirichlet=[],)
# spatial FESpace for level set (reference configuration)
#fes_lset_slice = H1(mesh, order=1, dirichlet=[])
# spatial FESpace for deformation field
fes_dfm_slice = H1(mesh, order=k_s, dim=mesh.dim)
# polynomial order in time for level set approximation
lset_order_time = 1
# integration order in time
time_order = 2
# time finite element (nodal!)
tfe = ScalarTimeFE(k_t) 
# space-time finite element space
st_fes = SpaceTimeFESpace(fes1,tfe, flags = {"dgjumps": True})

#Fitted heat equation example
tend = 0.5
delta_t = 1/64
tnew = 0

told = Parameter(0)
tref = ReferenceTimeVariable()
t = told + delta_t*tref

lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = lset_order_time,
                                         threshold=0.5, discontinuous_qn=True)

# radius of disk (the geometry)
r0 = 0.5

# position shift of the geometry in time
rho =  CoefficientFunction((1/(pi))*sin(2*pi*t))
# velocity of position shift
d_rho = CoefficientFunction(2*cos(2*pi*t))
#convection velocity:
w = CoefficientFunction((0,d_rho)) 

# level set
r = sqrt(x**2+(y-rho)**2)
levelset= r - r0

# solution and r.h.s.
Q = pi/r0   
u_exact = cos(Q*r) * sin(pi*t)
coeff_f = (Q/r * sin(Q*r) + (Q**2) * cos(Q*r)) * sin(pi*t) + pi * cos(Q*r) * cos(pi*t)


u0 = GridFunction(st_fes)

u0_ic = CreateTimeRestrictedGF(u0,0)

u,v = st_fes.TnT()


lset_p1 = lset_adap_st.lset_p1    

lset_top = CreateTimeRestrictedGF(lset_p1,1.0)
lset_bottom = CreateTimeRestrictedGF(lset_p1,0.0)

dfm = lset_adap_st.deform
dfm_top = CreateTimeRestrictedGF(dfm,1.0) # doesn't work for vector-valued functions yet...
#dfm_top = GridFunction(fes_dfm_slice)

t_old = 0
u0_ic.Set(u_exact)

h = specialcf.mesh_size

Draw(lset_top,mesh,"lset")
Draw(IfPos(-lset_top,u_exact,float('nan')),mesh,"u_exact")
Draw(IfPos(-lset_top,u0_ic,float('nan')),mesh,"u")
visoptions.deformation = 1

lset_neg = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}
lset_neg_bottom = { "levelset" : lset_bottom, "domain_type" : NEG, "subdivlvl" : 0}
lset_neg_top = { "levelset" : lset_top, "domain_type" : NEG, "subdivlvl" : 0}

def SpaceTimeNegBFI(form):
    return SymbolicBFI(levelset_domain = lset_neg, form = form, time_order=time_order, definedonelements = ci.GetElementsOfType(HASNEG))

ci = CutInfo(mesh,time_order=time_order)

    
hasneg_integrators_a = []
hasneg_integrators_f = []
patch_integrators_a = []

hasneg_integrators_a.append(SpaceTimeNegBFI(form = -u*dt(v)))
hasneg_integrators_a.append(SpaceTimeNegBFI(form = -delta_t*u*InnerProduct(w,grad(v))))
hasneg_integrators_a.append(SpaceTimeNegBFI(form = delta_t*grad(u)*grad(v)))
hasneg_integrators_a.append(SymbolicBFI(levelset_domain = lset_neg_top, form = fix_t(u,1)*fix_t(v,1)))
patch_integrators_a.append(SymbolicFacetPatchBFI(form = delta_t*0.01*h**(-2)*(u-u.Other())*(v-v.Other()),
                            skeleton=False,
                            time_order=time_order))
hasneg_integrators_f.append(SymbolicLFI(levelset_domain = lset_neg, form = delta_t*coeff_f*v, time_order=time_order)) 
hasneg_integrators_f.append(SymbolicLFI(levelset_domain = lset_neg_bottom,form = shifted_eval(u0_ic,back = dfm_top)*fix_t(v,0)))


a = BilinearForm(st_fes,check_unused=False,symmetric=False)
for integrator in hasneg_integrators_a + patch_integrators_a:
    a += integrator

f = LinearForm(st_fes)

for integrator in hasneg_integrators_f:
    f += integrator


while tend - t_old > delta_t/2:
    # update lset geometry to new time slab (also changes lset_p1 !)
    dfm = lset_adap_st.CalcDeformation(levelset,told,t_old,delta_t) 
    RestrictGFInTime(spacetime_gf=lset_p1,reference_time=0.0,space_gf=lset_bottom)
    RestrictGFInTime(spacetime_gf=lset_p1,reference_time=1.0,space_gf=lset_top)

    # update markers in (space-time) mesh
    ci.Update(lset_p1,time_order=time_order)

    # re-compute the facets for stabilization:
    ba_facets = GetFacetsWithNeighborTypes(mesh,a=ci.GetElementsOfType(HASNEG),
                                                b=ci.GetElementsOfType(IF))
    # re-evaluate the "active dofs" in the space time slab
    active_dofs = GetDofsOfElements(st_fes,ci.GetElementsOfType(HASNEG))

    # re-set definedonelements-markers according to new markings:
    for integrator in hasneg_integrators_a + hasneg_integrators_f:
        integrator.SetDefinedOnElements(ci.GetElementsOfType(HASNEG))
    for integrator in patch_integrators_a:
        integrator.SetDefinedOnElements(ba_facets)

    # update mesh deformation and assemble linear system
    mesh.SetDeformation(dfm)
    a.Assemble()
    f.Assemble()
    mesh.UnsetDeformation()

    # solve linear system
    u0.vec.data = a.mat.Inverse(active_dofs,"pardiso") * f.vec
       
    # evaluate upper trace of solution for
    #  * for error evaluation 
    #  * upwind-coupling to next time slab
    RestrictGFInTime(spacetime_gf=u0,reference_time=1.0,space_gf=u0_ic)   
    
    # update time variable (float and ParameterCL)
    t_old = t_old + delta_t
    told.Set(t_old)
    
    # compute error at final time
    mesh.SetDeformation(dfm_top)
    l2error = sqrt(Integrate(lset_neg_top,(u_exact-u0_ic)*(u_exact-u0_ic),mesh))
    mesh.UnsetDeformation()

    # print time and error
    print("\rt = {0:10}, l2error = {1:20}".format(t_old,l2error),end="")
    
    # Redraw:
    Redraw(blocking=True)

    # store deformation at top level of for next time step
    RestrictGFInTime(spacetime_gf=dfm,reference_time=1.0,space_gf=dfm_top)   

print("")       