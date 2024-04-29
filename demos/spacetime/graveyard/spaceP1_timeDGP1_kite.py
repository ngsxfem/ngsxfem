# unfitted Heat equation with Neumann b.c.
# solved with a P1-DG-in-time space-time discretization
from ngsolve import *
from netgen.geom2d import unit_square

from xfem import *
from math import pi
import sys

from xfem.lset_spacetime import *

use_sympy = True

if use_sympy:
    from kite_in_sympy import get_Laplace_u_str, get_dt_rho

ngsglobals.msg_level = 1

i = 3

if hasattr(sys, 'argv') and len(sys.argv) == 3 and sys.argv[1] == "i":
    i = int(sys.argv[2])
    print("Sys argv :", sys.argv)
    print("Loading manual val for i: ", i)
else:
    print("Loading default val for i: ", i)

square = SplineGeometry()
square.AddRectangle([-3.5,-1.5],[3.5,1.5])
ngmesh = square.GenerateMesh(maxh=(1/2)**i, quad_dominated=False)
mesh = Mesh (ngmesh)

coef_told = Parameter(0)
coef_delta_t = Parameter(0)
tref = ReferenceTimeVariable()
t = coef_told + coef_delta_t*tref

r0 = 1

rho =  (1 - 3*y**2)*t*sin(pi*t)
rhoL = lambda t: (1- 3*y**2)*t*sin(pi*t)
rho_str = "(1 - 3*y**2)*t*sin(pi*t)"

#convection velocity:
if use_sympy:
    d_rho = eval(get_dt_rho(rho_str))
else:
    rho = (1 - y**2)*t
    rhoL = lambda t: (1 - y**2)*t
    rho_str = "(1 - y**2)*t"
    d_rho = CoefficientFunction (1 - y**2)

w = CoefficientFunction((d_rho,0))

# level set
r = sqrt((x- rho)**2+y**2)
levelset= r - r0
levelsetL = lambda t: sqrt((x- rhoL(t))**2+y**2) - r0

Draw(levelsetL(0.), mesh,"lset", autoscale=False, min = 0, max = 0)

Q = pi/r0   
u_exact = cos(Q*r) * sin(pi*t)
u_exactL = lambda t: cos(Q*sqrt((x- rhoL(t))**2+y**2)) * sin(pi*t)
u_str = "cos(Q*r) * sin(pi*t)"

if use_sympy:
    Laplaceu = eval(get_Laplace_u_str(rho_str, u_str))
else:
    drdx = (x-rho)/r
    dr2dx2 = y**2/r**3

    drdy = y/r*(1 + 2*(x-rho)*t)
    dr2dy2 = (r - y*drdy)/r**2*(1+2*(x-rho)*t) + 4*t**2*y**2/r

    Laplaceu = - Q*sin(pi*t)*(cos(Q*r)*Q* ( drdx**2 + drdy**2) + sin(Q*r)*( dr2dx2 + dr2dy2 ))

coeff_f = -Laplaceu + pi * cos(Q*r) * cos(pi*t)
u_init = u_exact

# polynomial order in time
k_t = 1
# polynomial order in space
k_s = 1
# spatial FESpace for solution
fes1 = H1(mesh, order=k_s)
# polynomial order in time for level set approximation
lset_order_time = 1
# integration order in time
time_order = 2
# time finite element (nodal!)
tfe = ScalarTimeFE(k_t)
# space-time finite element space
st_fes = SpaceTimeFESpace(fes1,tfe, flags = {"dgjumps": True})

#Unfitted heat equation example
tend = 1
delta_t = tend/(2**(i+2))
coef_delta_t.Set(delta_t)
tnew = 0
told = 0

lset_p1 = GridFunction(st_fes)

SpaceTimeInterpolateToP1(levelset,tref,lset_p1)

lset_top = CreateTimeRestrictedGF(lset_p1,1.0)
lset_bottom = CreateTimeRestrictedGF(lset_p1,0.0)

gfu = GridFunction(st_fes)

u_last = CreateTimeRestrictedGF(gfu,0)
u_last.Set(u_exactL(0.))

u,v = st_fes.TnT()

h = specialcf.mesh_size

lset_neg = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}
lset_neg_bottom = { "levelset" : lset_bottom, "domain_type" : NEG, "subdivlvl" : 0}
lset_neg_top = { "levelset" : lset_top, "domain_type" : NEG, "subdivlvl" : 0}

def SpaceTimeNegBFI(form):
    return SymbolicBFI(levelset_domain = lset_neg, form = form, time_order=time_order)

ci = CutInfo(mesh,time_order=time_order)

#n_F = specialcf.normal(mesh.dim)

hasneg_integrators_a = []
hasneg_integrators_f = []
patch_integrators_a = []
hasneg_integrators_a.append(SpaceTimeNegBFI(form = delta_t*grad(u)*grad(v)))
hasneg_integrators_a.append(SymbolicBFI(levelset_domain = lset_neg_top, form = fix_t(u,1)*fix_t(v,1)))
hasneg_integrators_a.append(SpaceTimeNegBFI(form = -u*dt(v)))
hasneg_integrators_a.append(SpaceTimeNegBFI(form = -delta_t*u*InnerProduct(w,grad(v))))
patch_integrators_a.append(SymbolicFacetPatchBFI(form = delta_t*1.05*h**(-2)*(u-u.Other())*(v-v.Other()),
                                                 skeleton=False, time_order=time_order))
#patch_integrators_a.append(SymbolicFacetPatchBFI( InnerProduct(grad(u) - grad(u.Other()), n_F) * InnerProduct(grad(v) - grad(v.Other()), n_F), skeleton=True, time_order=time_order))
hasneg_integrators_f.append(SymbolicLFI(levelset_domain = lset_neg, form = delta_t*coeff_f*v, time_order=time_order)) 
hasneg_integrators_f.append(SymbolicLFI(levelset_domain = lset_neg_bottom,form = u_last*fix_t(v,0)))


a = BilinearForm(st_fes,check_unused=False,symmetric=False)
for integrator in hasneg_integrators_a + patch_integrators_a:
    a += integrator

f = LinearForm(st_fes)

for integrator in hasneg_integrators_f:
    f += integrator

while tend - told > delta_t/2:
    SpaceTimeInterpolateToP1(levelset,tref,lset_p1)
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

    # assemble linear system
    a.Assemble()
    f.Assemble()

    # solve linear system
    inv = a.mat.Inverse(active_dofs,inverse="")
    gfu.vec.data =  inv * f.vec
       
    # evaluate upper trace of solution for
    #  * for error evaluation 
    #  * upwind-coupling to next time slab
    RestrictGFInTime(spacetime_gf=gfu,reference_time=1.0,space_gf=u_last)   

    # update time variable (float and ParameterCL)
    told = told + delta_t
    coef_told.Set(told)
    
    Draw(levelsetL(told), mesh,"lset", autoscale=False, min = 0, max = 0)
    #Draw(IfPos(-levelsetL(told), u_exactL(told) ,float('nan')),mesh,"u", autoscale = False, min = -1, max = 1)
    #Draw(IfPos(-levelsetL(told), u_exactL(told) ,float('nan')),mesh,"u", autoscale = False, min = -1, max = 1)
    #Draw(IfPos(-levelsetL(told), u_exactL(told)-u_last, float('nan')),mesh,"error")
    
    # compute error at end of time slab
    l2error = sqrt(Integrate(lset_neg_top,(u_exactL(told) -u_last)**2,mesh))
    # print time and error
    print("t = {0:10}, l2error = {1:20}".format(told,l2error),end="\n")
