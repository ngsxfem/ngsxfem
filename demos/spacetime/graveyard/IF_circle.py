from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from numpy import pi

from xfem.lset_spacetime import *

# geometry
# Geometry 
square = SplineGeometry()
square.AddRectangle([0,0],[2,2],bc=1)
maxh = 0.1
mesh = Mesh (square.GenerateMesh(maxh=maxh, quad_dominated=False))

n_ref = 0
for i in range(n_ref):
    mesh.Refine()

h = specialcf.mesh_size
n_outer = specialcf.normal(mesh.dim)



# FE Spaces
k_t = 2
k_s = 2
fes1 = H1(mesh, order=k_s,dirichlet=[1])
fes_lset_slice = H1(mesh, order=1, dirichlet=[])
fes_dfm_slice = H1(mesh, order=k_s, dim=mesh.dim)

lset_order_time = 2
time_order = 4
tfe = ScalarTimeFE(k_t) 

st_fes = SpaceTimeFESpace(fes1,tfe,dirichlet=[1])
VhG_ic = FESpace([fes1,fes1])
VhG = FESpace([st_fes,st_fes])


tend = 0.5
delta_t = 1/16
tnew = 0

told = Parameter(0)
tref = ReferenceTimeVariable()
t = told + delta_t*tref

lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = lset_order_time,
                                         threshold=0.1, discontinuous_qn=True)

# Problem data
rho =  CoefficientFunction((1/(4*pi))*sin(2*pi*t))
d_rho = CoefficientFunction(0.5*cos(2*pi*t))

#rho =  CoefficientFunction(0)
#d_rho = CoefficientFunction(0)

alpha_neg = 10.0
alpha_pos = 20.0
alpha = [alpha_neg,alpha_pos]
beta_neg = 2.0
beta_pos = 1.0
beta = [beta_neg,beta_pos]

r0 = 1/3
w = CoefficientFunction((d_rho,0)) # convection
levelset= CoefficientFunction(sqrt((x-0.5-rho)*(x-0.5-rho)+(y-1)*(y-1)) -r0)
a_const = 1.1569
b_const = -8.1621
k_const = 1

r_lset = CoefficientFunction(sqrt((x-0.5-rho)*(x-0.5-rho) + (y-1)*(y-1)))
f_pos = CoefficientFunction(   pi*cos(pi*t)*cos(pi*r_lset) 
                             + alpha_pos*sin(pi*t)*pi*pi*cos(pi*r_lset)
                             + alpha_pos*pi*sin(pi*t)*sin(pi*r_lset)/r_lset )
f_neg = ( k_const*pi*cos(k_const*pi*t)*(a_const+b_const*((x-0.5-rho)*(x-0.5-rho)+(y-1)*(y-1))) 
         - alpha_neg*sin(k_const*pi*t)*2*b_const
         - alpha_neg*sin(k_const*pi*t)*2*b_const   )
coeff_f = [f_neg,f_pos]


# for monitoring the error
exact = [sin(k_const*pi*t)*(a_const+b_const*((x-0.5-rho)*(x-0.5-rho)+(y-1)*(y-1)) ),
         sin(k_const*pi*t)*cos(pi*sqrt((x-0.5-rho)*(x-0.5-rho)+(y-1)*(y-1)))]

u0 = GridFunction(VhG)
u0_ic = GridFunction(VhG_ic)
#Draw(u0_ic,mesh,"u")

# Definition of Trial/Test-functions
u = VhG.TrialFunction()
v = VhG.TestFunction()

gradu = [grad(ui) for ui in u]
gradv = [grad(vi) for vi in v]


lset_top = GridFunction(fes_lset_slice)
lset_bottom = GridFunction(fes_lset_slice)
dfm_top = GridFunction(fes_dfm_slice)

t_old = 0
coef_u0 = [CoefficientFunction(0),CoefficientFunction(0)]
u0_ic.components[0].Set(coef_u0[0]) # set initial condition
u0_ic.components[1].Set(coef_u0[1]) # set initial condition

stab = 10*(alpha[0]+alpha[1])*(k_s+1)*(k_s+1)/h
stab_outer = 10*(k_s+1)*(k_s+1)/h


while tend - t_old > delta_t/2:
    
    dfm = lset_adap_st.CalcDeformation(levelset,told,t_old,delta_t,calc_kappa = True) 
    dfm_top.vec[:] = dfm.vec[lset_order_time*lset_adap_st.ndof_node : (lset_order_time+1)*lset_adap_st.ndof_node]
    lset_p1 = lset_adap_st.lset_p1  
    active_dofs = BitArray(VhG.FreeDofs())
   
    hasneg_spacetime = lset_adap_st.hasneg_spacetime
    hasneg_spacetime |= lset_adap_st.hasif_spacetime
    haspos_spacetime = lset_adap_st.haspos_spacetime
    haspos_spacetime |= lset_adap_st.hasif_spacetime
    
    active_dofs &= CompoundBitArray([GetDofsOfElements(st_fes,hasneg_spacetime),GetDofsOfElements(st_fes,haspos_spacetime)])
    #print(active_dofs) 
    
    lset_neg = { "levelset" : lset_p1, "domain_type" : NEG, "subdivlvl" : 0}
    lset_pos = { "levelset" : lset_p1, "domain_type" : POS, "subdivlvl" : 0}
    lset_if = { "levelset" : lset_p1, "domain_type" : IF, "subdivlvl" : 0}
    # Put this into LevelSetMeshAdaptation_Spacetime later on
    lset_bottom.vec[:] = lset_p1.vec[0:fes_lset_slice.ndof]
    lset_neg_bottom = { "levelset" : lset_bottom, "domain_type" : NEG, "subdivlvl" : 0}
    lset_pos_bottom = { "levelset" : lset_bottom, "domain_type" : POS, "subdivlvl" : 0}
    lset_top.vec[:] = lset_p1.vec[lset_order_time*fes_lset_slice.ndof : (lset_order_time+1)*fes_lset_slice.ndof]
    lset_neg_top = { "levelset" : lset_top, "domain_type" : NEG, "subdivlvl" : 0}
    lset_pos_top = { "levelset" : lset_top, "domain_type" : POS, "subdivlvl" : 0}
    
    n_lset_p1 = 1.0/grad(lset_adap_st.lset_p1).Norm() * grad(lset_adap_st.lset_p1)
    
    
    kappas = {"interp" : (lset_adap_st.kappa,1-lset_adap_st.kappa), 
              "st_CutInfo" : (CutRatioGF(lset_adap_st.ci),1-CutRatioGF(lset_adap_st.ci)) }
    kappa = kappas["st_CutInfo"]
    
#    time_quad = lset_adap_st.v_ho_st.TimeFE_nodes().NumPy()
#    times = [t_old + delta_t * xi for xi in time_quad]
#    for tau in times:  
#        lset_adap_st.v_kappa.SetTime(tau)
#        Draw(kappa[0],mesh,"kappa")
#        input("Kappa at t = {0}".format(tau))
#    lset_adap_st.v_kappa.SetOverrideTime(False)
         
    average_flux_u = sum( [-kappa[i]*alpha[i]*gradu[i]*n_lset_p1 for i in [0,1]] )
    average_flux_v = sum( [-kappa[i]*alpha[i]*gradv[i]*n_lset_p1 for i in [0,1]] )
    betajump_u = beta[0]*u[0]-beta[1]*u[1]
    betajump_v = beta[0]*v[0]-beta[1]*v[1]
    
     
    a = BilinearForm(VhG,symmetric=False)

    a += SymbolicBFI(levelset_domain = lset_neg, form =  -u[0]*dt(v[0],0)*beta[0]
                                                      + delta_t*alpha[0]*beta[0]*gradu[0]*gradv[0]
                                                      - delta_t*beta[0]*u[0]*gradv[0]*w, time_order=time_order)    
                                                      
    a += SymbolicBFI(levelset_domain = lset_pos, form = -u[1]*dt(v[1],1)*beta[1]  
                                                      + delta_t*alpha[1]*beta[1]*gradu[1]*gradv[1]
                                                      - delta_t*beta[1]*u[1]*gradv[1]*w, time_order=time_order)
 
    
    a += SymbolicBFI(levelset_domain = lset_neg_top, form = beta[0]*fix_t(u[0],1,0)*fix_t(v[0],1,0) )
    a += SymbolicBFI(levelset_domain = lset_pos_top, form = beta[1]*fix_t(u[1],1,1)*fix_t(v[1],1,1) )
    
    

    a += SymbolicBFI(levelset_domain = lset_if, form = delta_t*( average_flux_u*betajump_v   
                                                              + betajump_u*average_flux_v 
                                                            + stab*betajump_u*betajump_v), time_order=time_order )

                                                               
                
    f = LinearForm(VhG)
    f += SymbolicLFI(levelset_domain = lset_neg, form = delta_t*beta[0]*coeff_f[0]*v[0], time_order=time_order)
    f += SymbolicLFI(levelset_domain = lset_pos, form = delta_t*beta[1]*coeff_f[1]*v[1], time_order=time_order)
    f += SymbolicLFI(levelset_domain = lset_neg_bottom,form = beta[0]*u0_ic.components[0]*fix_t(v[0],0,0) )
    f += SymbolicLFI(levelset_domain = lset_pos_bottom,form = beta[1]*u0_ic.components[1]*fix_t(v[1],0,1) )
     
    
    mesh.SetDeformation(dfm)
    a.Assemble()
    f.Assemble()

    mesh.UnsetDeformation()
    
    # Setting boundary values (temporary solution)
    u0.vec[:] = 0
   
    times = [t_old + delta_t * xi for xi in st_fes.TimeFE_nodes()]
    u_bnd = GridFunction(fes1)
    for i,ti in enumerate(times):
        rho_ti = (1/(4*pi))*sin(2*pi*ti)
        print(ti)
        u_bnd.Set(sin(k_const*pi*ti)*cos(pi*sqrt((x-0.5-rho_ti)*(x-0.5-rho_ti)+(y-1)*(y-1))),BND)
        u0.components[1].vec[i*fes1.ndof : (i+1)*fes1.ndof] = u_bnd.vec[:] 
    
    #u_bnd.vec[:] = u0.components[1].vec[0 : fes1.ndof]   
    #Draw(sqrt((u_bnd-exact[1])*(u_bnd-exact[1])),mesh,"errorB")
    #input("Stop")
    
           
    res = f.vec.CreateVector()
    res.data = f.vec - a.mat * u0.vec
    u0.vec.data += a.mat.Inverse(active_dofs,"") * res
           
    
       
    # exploiting the nodal property of the time fe:b_const
    for i in range(2): 
        u0_ic.components[i].vec[:] = u0.components[i].vec[k_t*fes1.ndof : (k_t+1)*fes1.ndof]
   
    
    
    t_old = t_old + delta_t
    told.Set(t_old)
    #u0_ic.components[1].Set(exact[1])
    
    mesh.SetDeformation(dfm_top)
    err_sqr_coefs = [ (u0_ic.components[i] - exact[i])*(u0_ic.components[i] - exact[i]) for i in [0,1] ]
    l2error = sqrt( Integrate( levelset_domain=lset_neg_top, cf=err_sqr_coefs[0], mesh=mesh,
                           order=2*k_s, heapsize=1000000)
                + Integrate(levelset_domain=lset_pos_top, cf=err_sqr_coefs[1], mesh=mesh,
                            order=2*k_s, heapsize=1000000)) 
    mesh.UnsetDeformation()
    print("t = {0}, l2error = {1}".format(t_old,l2error))
    
    Draw(IfPos(lset_top,u0_ic.components[1],u0_ic.components[0]),mesh,"u")
    #Draw(IfPos(lset_top,sqrt(err_sqr_coefs[1]),sqrt(err_sqr_coefs[0])),mesh,"error")
    #Draw(IfPos(lset_top,exact[1],exact[0]),mesh,"exact")
    #input("Continue")
    #Redraw(blocking=True)
       
