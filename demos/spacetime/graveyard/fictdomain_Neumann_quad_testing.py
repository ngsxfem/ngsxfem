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
ngmesh = square.GenerateMesh(maxh=10, quad_dominated=True)
mesh = Mesh (ngmesh)

for i in range(6):
    mesh.Refine()

fes1 = H1(mesh, order=1, dirichlet=[])
fes_lset_slice = H1(mesh, order=1, dirichlet=[])
k_t = 1
k_s = 1
lset_order_time = 2
tfe = ScalarTimeFE(k_t) 

st_fes = SpaceTimeFESpace(fes1,tfe)
visoptions.deformation = 1

#Fitted heat equation example
tend = 0.5
delta_t = 1/32
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
Draw(u0_ic,mesh,"u")
u = st_fes.TrialFunction()
v = st_fes.TestFunction()

lset_top = GridFunction(fes_lset_slice)
lset_bottom = GridFunction(fes_lset_slice)

t_old = 0
u0_ic.Set(u_exact)


while tend - t_old > delta_t/2:
    
    dfm = lset_adap_st.CalcDeformation(levelset,told,t_old,delta_t) 
    lset_p1 = lset_adap_st.lset_p1    
    hasneg_spacetime = lset_adap_st.hasneg_spacetime
    active_dofs = GetDofsOfElements(st_fes,hasneg_spacetime)
    #print(active_dofs)
    
    lset_neg = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}
    # Put this into LevelSetMeshAdaptation_Spacetime later on
    lset_bottom.vec[:] = lset_p1.vec[0:fes_lset_slice.ndof]
    lset_neg_bottom = { "levelset" : lset_bottom, "domain_type" : NEG, "subdivlvl" : 0}
    lset_top.vec[:] = lset_p1.vec[lset_order_time*fes_lset_slice.ndof : (lset_order_time+1)*fes_lset_slice.ndof]
    lset_neg_top = { "levelset" : lset_top, "domain_type" : NEG, "subdivlvl" : 0}
    
    a = BilinearForm(st_fes,symmetric=False)
    a += SymbolicBFI(levelset_domain = lset_neg, form = -u*dt(v), time_order=2)
    a += SymbolicBFI(levelset_domain = lset_neg, form = delta_t*grad(u)*grad(v), time_order=2)
    a += SymbolicBFI(levelset_domain = lset_neg, form = -delta_t*u*w*grad(v), time_order=2)
    a += SymbolicBFI(levelset_domain = lset_neg_top, form = fix_t(u,1)*fix_t(v,1) )
    
    #Draw(lset_bottom,mesh,"bottom") 
    #Draw(lset_top,mesh,"top") 
    #input("lvlsets")
    a.Assemble()
              
    f = LinearForm(st_fes)
    f += SymbolicLFI(levelset_domain = lset_neg, form = delta_t*coeff_f*v, time_order=2)
    f += SymbolicLFI(levelset_domain = lset_neg_bottom,form = u0_ic*fix_t(v,0) )
    f.Assemble()
    u0.vec.data = a.mat.Inverse(active_dofs,"") * f.vec
       
    # exploiting the nodal property of the time fe:
    u0_ic.vec[:] = u0.vec[fes1.ndof : 2*fes1.ndof]
    
    t_old = t_old + delta_t
    told.Set(t_old)
    
    l2error = sqrt(Integrate(lset_neg_top,(u_exact-u0_ic)*(u_exact-u0_ic),mesh))
    print("t = {0}, l2error = {1}".format(t_old,l2error))
   
    Redraw(blocking=True)
    #Draw(sqrt((u_exact-u0_ic)*(u_exact-u0_ic)),mesh,"error") 
       
