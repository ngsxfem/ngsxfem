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
k_s = 1
# spatial FESpace for solution
fes1 = H1(mesh, order=k_s, dirichlet=[],)
# spatial FESpace for level set (reference configuration)
fes_lset_slice = H1(mesh, order=1, dirichlet=[])
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
delta_t = 1/32
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
u0_ic = GridFunction(fes1)

u,v = st_fes.TnT()

lset_top = GridFunction(fes_lset_slice)
lset_bottom = GridFunction(fes_lset_slice)
dfm_top = GridFunction(fes_dfm_slice)

t_old = 0
u0_ic.Set(u_exact)

h = specialcf.mesh_size

Draw(lset_top,mesh,"lset")
Draw(IfPos(-lset_top,u_exact,float('nan')),mesh,"u_exact")
Draw(IfPos(-lset_top,u0_ic,float('nan')),mesh,"u")
visoptions.deformation = 1

while tend - t_old > delta_t/2:
    
    dfm = lset_adap_st.CalcDeformation(levelset,told,t_old,delta_t) 
    dfm_top.vec[:] = dfm.vec[lset_order_time*lset_adap_st.ndof_node : (lset_order_time+1)*lset_adap_st.ndof_node]
    lset_p1 = lset_adap_st.lset_p1    

    ci = CutInfo(mesh,lset_p1,time_order=time_order)

    ba_facets = GetFacetsWithNeighborTypes(mesh,a=ci.GetElementsOfType(HASNEG),
                                                b=ci.GetElementsOfType(IF))

    active_dofs = GetDofsOfElements(st_fes,ci.GetElementsOfType(HASNEG))
    
    lset_neg = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}

    # Put this into LevelSetMeshAdaptation_Spacetime later on
    lset_bottom.vec[:] = lset_p1.vec[0:fes_lset_slice.ndof]
    lset_neg_bottom = { "levelset" : lset_bottom, "domain_type" : NEG, "subdivlvl" : 0}
    lset_top.vec[:] = lset_p1.vec[lset_order_time*fes_lset_slice.ndof : (lset_order_time+1)*fes_lset_slice.ndof]
    lset_neg_top = { "levelset" : lset_top, "domain_type" : NEG, "subdivlvl" : 0}
    
    def SpaceTimeNegBFI(form):
        return SymbolicBFI(levelset_domain = lset_neg, form = form, time_order=time_order, definedonelements = ci.GetElementsOfType(HASNEG))

    a = BilinearForm(st_fes,check_unused=False,symmetric=False)
    a += SpaceTimeNegBFI(form = -u*dt(v))
    a += SpaceTimeNegBFI(form = delta_t*grad(u)*grad(v))
    a += SpaceTimeNegBFI(form = -delta_t*u*InnerProduct(w,grad(v)))
    a += SymbolicBFI(levelset_domain = lset_neg_top, form = fix_t(u,1)*fix_t(v,1), definedonelements = ci.GetElementsOfType(HASNEG) )
    a += SymbolicFacetPatchBFI(form = delta_t*0.01*h**(-2)*(u-u.Other())*(v-v.Other()),
                               skeleton=False,
                               definedonelements=ba_facets,
                               time_order=time_order)

    f = LinearForm(st_fes)
    f += SymbolicLFI(levelset_domain = lset_neg, form = delta_t*coeff_f*v, time_order=time_order, definedonelements = ci.GetElementsOfType(HASNEG))
    f += SymbolicLFI(levelset_domain = lset_neg_bottom,form = u0_ic*fix_t(v,0),definedonelements = ci.GetElementsOfType(HASNEG) )
    
    mesh.SetDeformation(dfm)
    a.Assemble()
    f.Assemble()
    mesh.UnsetDeformation()

    u0.vec.data = a.mat.Inverse(active_dofs,"pardiso") * f.vec
       
    u0_ic.vec[:] = u0.vec[k_t*fes1.ndof : (k_t+1)*fes1.ndof]
    
    t_old = t_old + delta_t
    told.Set(t_old)
    
    mesh.SetDeformation(dfm_top)
    l2error = sqrt(Integrate(lset_neg_top,(u_exact-u0_ic)*(u_exact-u0_ic),mesh))
    mesh.UnsetDeformation()
    print("\rt = {0:10}, l2error = {1:20}".format(t_old,l2error),end="")
    
    Redraw(blocking=True)
print("")       