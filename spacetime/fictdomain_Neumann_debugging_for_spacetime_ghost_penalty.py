# unfitted Heat equation with Neumann b.c.
from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from numpy import pi

from xfem.lset_spacetime import *

square = SplineGeometry()
square.AddRectangle([-1,-0.75],[1,1.5],bc=1)
ngmesh = square.GenerateMesh(maxh=0.1, quad_dominated=False)
mesh = Mesh (ngmesh)

n_ref = 0
for i in range(n_ref):
    mesh.Refine()

k_t = 1
k_s = 2
fes1 = H1(mesh, order=k_s, dirichlet=[])
fes_lset_slice = H1(mesh, order=1, dirichlet=[])
fes_dfm_slice = H1(mesh, order=k_s, dim=mesh.dim)

lset_order_time = 2
time_order = 2
tfe = ScalarTimeFE(k_t) 

st_fes = SpaceTimeFESpace(fes1,tfe, flags = {"dgjumps": True})
visoptions.deformation = 1

#Fitted heat equation example
tend = 0.5
delta_t = 1/64
tnew = 0

told = Parameter(0)
tref = ReferenceTimeVariable()
t = told + delta_t*tref

lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = lset_order_time,
                                         threshold=0.5, discontinuous_qn=True)

r0 = 0.5
tmp1 = 0.5*pi/r0                                         
rho =  CoefficientFunction((1/(pi))*sin(2*pi*t))
d_rho = CoefficientFunction(2*cos(2*pi*t))
w = CoefficientFunction((0,d_rho)) # convection
levelset= CoefficientFunction(sqrt(x*x+(y-rho)*(y-rho)) -r0)
coeff_f = CoefficientFunction(  -(pi/r0)*tmp1*( sin(tmp1*sqrt(x*x+(y-rho)*(y-rho)))*sin(tmp1*sqrt(x*x+(y-rho)*(y-rho))) 
                                - cos(tmp1*sqrt(x*x+(y-rho)*(y-rho)))*cos(tmp1*sqrt(x*x+(y-rho)*(y-rho))) )  
                                + (pi/r0)*cos(tmp1*sqrt(x*x+(y-rho)*(y-rho)))*sin(tmp1*sqrt(x*x+(y-rho)*(y-rho)))*(1/sqrt(x*x+(y-rho)*(y-rho))) )

# for monitoring the error
u_exact = CoefficientFunction(cos(tmp1*sqrt(x*x+(y-rho)*(y-rho)))*cos(tmp1*sqrt(x*x+(y-rho)*(y-rho))) )


u0 = GridFunction(st_fes)
u0_ic = GridFunction(fes1)
#Draw(u0_ic,mesh,"u")
u = st_fes.TrialFunction()
v = st_fes.TestFunction()

lset_top = GridFunction(fes_lset_slice)
lset_bottom = GridFunction(fes_lset_slice)
dfm_top = GridFunction(fes_dfm_slice)

t_old = 0
u0_ic.Set(u_exact)


# normal derivative jumps:
def dnjump(u,order,comp = -1):
    if order%2==0:
        return dn(u,order,comp) - dn(u.Other(),order,comp)
    else:
        return dn(u,order,comp) + dn(u.Other(),order,comp)


while tend - t_old > delta_t/2:
    
    dfm = lset_adap_st.CalcDeformation(levelset,told,t_old,delta_t) 
    dfm_top.vec[:] = dfm.vec[lset_order_time*lset_adap_st.ndof_node : (lset_order_time+1)*lset_adap_st.ndof_node]
    lset_p1 = lset_adap_st.lset_p1    

    hasneg_spacetime = lset_adap_st.hasneg_spacetime
    hasneg_spacetime |= lset_adap_st.hasif_spacetime
    active_dofs = GetDofsOfElements(st_fes,hasneg_spacetime)
    #print(active_dofs)
    
    lset_neg = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}
    # Put this into LevelSetMeshAdaptation_Spacetime later on
    lset_bottom.vec[:] = lset_p1.vec[0:fes_lset_slice.ndof]
    lset_neg_bottom = { "levelset" : lset_bottom, "domain_type" : NEG, "subdivlvl" : 0}
    lset_top.vec[:] = lset_p1.vec[lset_order_time*fes_lset_slice.ndof : (lset_order_time+1)*fes_lset_slice.ndof]
    lset_neg_top = { "levelset" : lset_top, "domain_type" : NEG, "subdivlvl" : 0}
    
    a = BilinearForm(st_fes,symmetric=False)
    a += SymbolicBFI(levelset_domain = lset_neg, form = -u*dt(v), time_order=time_order)
    a += SymbolicBFI(levelset_domain = lset_neg, form = delta_t*grad(u)*grad(v), time_order=time_order)
    a += SymbolicBFI(levelset_domain = lset_neg, form = -delta_t*u*w*grad(v), time_order=time_order)
    a += SymbolicBFI(levelset_domain = lset_neg_top, form = fix_t(u,1)*fix_t(v,1) )

    ba_facets = GetFacetsWithNeighborTypes(mesh,
                                           a=hasneg_spacetime,
                                           b=lset_adap_st.hasif_spacetime)

    # a += SymbolicBFI(form = 0.1 * specialcf.mesh_size * dnjump(fix_t(u,0),1) * dnjump(fix_t(v,0),1), VOL_or_BND = VOL, skeleton=True)
    # a += SymbolicBFI(form = 0.1 * specialcf.mesh_size * dnjump(fix_t(u,1),1) * dnjump(fix_t(v,1),1), VOL_or_BND = VOL, skeleton=True)

    # a += SymbolicBFI(levelset_domain = lset_neg_top, form = 1e-6 * fix_t(u,1) * fix_t(v,1))
    # a += SymbolicBFI(levelset_domain = lset_neg_bottom, form = 1e-6 * fix_t(u,0) * fix_t(v,0))
    
    # a += SymbolicBFI(form = 1e-6 * fix_t(u,0) * fix_t(v,0))
    # a += SymbolicBFI(form = 1e-6 * fix_t(u,1) * fix_t(v,1))
    
    
    a += SymbolicFacetPatchBFI(form = delta_t * specialcf.mesh_size * dnjump(u,1) * dnjump(v,1),
                               time_order=time_order,
                               definedonelements=ba_facets)
    a += SymbolicFacetPatchBFI(form = delta_t * specialcf.mesh_size *specialcf.mesh_size *specialcf.mesh_size * dnjump(u,2) * dnjump(v,2),
                               time_order=time_order,
                               definedonelements=ba_facets)
    # a += SymbolicFacetPatchBFI(form = delta_t * specialcf.mesh_size *  (grad(u)-grad(u.Other())) * (grad(v)-grad(v.Other())),
    #                            time_order=time_order,
    #                            definedonelements=ba_facets)
    # a += SymbolicFacetPatchBFI( (u-u.Other()) * (v-v.Other()),
    #                            time_order=time_order,
    #                            definedonelements=ba_facets)
    
    #Draw(lset_bottom,mesh,"bottom") 
    #Draw(lset_top,mesh,"top") 
    #input("lvlsets")
                
    f = LinearForm(st_fes)
    f += SymbolicLFI(levelset_domain = lset_neg, form = delta_t*coeff_f*v, time_order=time_order)
    f += SymbolicLFI(levelset_domain = lset_neg_bottom,form = u0_ic*fix_t(v,0) )
    
    mesh.SetDeformation(dfm)
    a.Assemble()
    f.Assemble()
    mesh.UnsetDeformation()
    u0.vec.data = a.mat.Inverse(active_dofs,"pardiso") * f.vec
    # u0.vec.data = a.mat.Inverse(st_fes.FreeDofs(),"umfpack") * f.vec
    # u0.vec.data = a.mat.Inverse(st_fes.FreeDofs(),"pardiso") * f.vec
       
    # exploiting the nodal property of the time fe:
    u0_ic.vec[:] = u0.vec[k_t*fes1.ndof : (k_t+1)*fes1.ndof]
    
    t_old = t_old + delta_t
    told.Set(t_old)
    
    mesh.SetDeformation(dfm_top)
    l2error = sqrt(Integrate(lset_neg_top,(u_exact-u0_ic)*(u_exact-u0_ic),mesh))
    mesh.UnsetDeformation()
    print("t = {0}, l2error = {1}".format(t_old,l2error))
    
    # Draw(IfPos(-lset_top,u0_ic,float('nan')),mesh,"u")
    # Draw(u0_ic,mesh,"u2")
    Draw(IfPos(BitArrayCF(hasneg_spacetime)-0.5,u0_ic,float('nan')),mesh,"u")
    #Redraw(blocking=True)
    #Draw(sqrt((u_exact-u0_ic)*(u_exact-u0_ic)),mesh,"error") 
       
