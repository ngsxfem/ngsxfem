"""
In this example we solve an *unfitted* Stokes interface problem with a high order isoparametric
unfitteddiscretization method. 
This file is based on the discretizations in cutfem.py, nxfem.py and nxfem_higher_order.py. 

    domain: 
    -------
    The domain is [-1,1]^2 while the interface is described by a level set function ( circile with
    radius R=2/3).

    PDE problem:
    ------------
    domain equations for velocity (u1,u2) and pressure p:
     - alpha_i (u1_xx + u1_yy) + p_x = rho_i g1 in subdomain i, i=1,2,
     - alpha_i (u2_xx + u2_yy) + p_y = rho_i g2 in subdomain i, i=1,2,
                         u1_x + u2_y = 0        in subdomain i, i=1,2,
    interface conditions:
                                [u] =    0 on interface (continuity across the interface    ),
             [-alpha · du/dn + p·n] =    f on interface (conservation of the (momentum) flux),
                                 u  =  u_D on domain boundary.

    The r.h.s. term (g1,g2) corresponds to gravity, the term f is surface tension force (here f =
    kappa · n where kappa is the mean curvature (1/R)). The Dirichlet data is chosen according to a
    manufactured solution introduced in [1]
    which allows us to measure errors after the computation of a discrete solution.
    The coefficients alpha are domain-wise constants which are different in the two subdomains.
    

    discretization:
    ---------------
    Finite element space:
    As in nxfem.py but for every velocity component and the pressure

    Variational formulation:
    We use a Nitsche formulation which involves averages of the fluxes and jumps of the solution
    across the interface [2]. For the average we use the Hansbo-choice [3] where the average is
    adjusted to the local cut configuration in order to ensure stability of the resulting viscosity
    formulation. 
    To ensure inf-sup for the velocity-pressure space we add a ghost penalty stabilization on the
    pressure space, cf. [2]. 

    Surface tension:  
    In this example we prescribe the surface tension analytically, i.e. we circumvent approximating
    the mean curvature from the level set function. 

    implementational aspects:
    ---------------
    Geometry approximation:
    As in nxfem_higher_order.py

    Ghost penalty stabilization:
    The edge-based stabilizations require different couplings (compared to standard
    discretizations). To add these to the sparsity pattern we have to add the "dgjumps" flags which
    prepares the sparse matrix for the corresponding needed couplings.


    linear systems:
    ---------------
    A (sparse) direct solver is applied to solve the arising linear systems.

    literature:
    -----------
    [1]: M.Kirchhart, S.Groß, A.Reusken, Analysis of an XFEM discretization for Stokes interface
    problems, SISC, 2016
    [2]: P.Lederer, C.-M.Pfeiler, C.Wintersteiger, C.Lehrenfeld, Higher order unfitted fem for
    stokes interface problems. PAMM, 2016 
    [3]: A.Hansbo, P.Hansbo, An unfitted finite element method, based on Nitsche's method, for
    elliptic interface problems, Comp. Meth. Appl. Mech. Eng., 2002

"""

from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters
from ngsolve import *
from xfem import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *
from math import pi
# 2D: circle configuration
square = SplineGeometry()
square.AddRectangle([-1,-1],[1,1],bcs=[1,2,3,4])
mesh = Mesh (square.GenerateMesh(maxh=4, quad_dominated=False))
for i in range(3):
    mesh.Refine()

order = 2 #=k in P(k)P(k-1)-Taylor-Hood pair
    
mu1 = 1.0
mu2 = 10.0
mu = [mu1,mu2]

R = 2.0/3.0
aneg = 1.0/mu1
apos = 1.0/mu2 + (1.0/mu1 - 1.0/mu2)*exp(x*x+y*y-R*R)
gammaf = 0.5 # surface tension = pressure jump
#q = gammaf - pi*R*R/4.0*gammaf

functions = {
           "Levelset" : sqrt(x*x+y*y)-R,
           "SourceInnerX" : (exp(-1*(x*x+y*y))*((-8*y)+(4*x*x*y)+(4*y*y*y))+3*x*x),
           "SourceOuterX" : (exp(-1*(x*x+y*y))*((-8*y)+(4*x*x*y)+(4*y*y*y))+3*x*x),
           "SourceInnerY" : (exp(-1*(x*x+y*y))*((-4*x*x*x)+(8*x)-(4*x*y*y))),
           "SourceOuterY" : (exp(-1*(x*x+y*y))*((-4*x*x*x)+(8*x)-(4*x*y*y))),
           "SolutionOuterVelX" : (apos*exp(-1.0 *( x * x + y * y)) * -1.0 * y),
           "SolutionOuterVelY" : (apos*exp(-1.0 *( x * x + y * y)) * x),
           "SolutionInnerVelX" : (aneg*exp(-1.0 *( x * x + y * y)) * -1.0 * y),
           "SolutionInnerVelY" : (aneg*exp(-1.0 *( x * x + y * y)) * x),
           "SolutionOuterVelX_DX" : (2.0*x*y/mu2*exp(-1.0 *( x * x + y * y))),
           "SolutionOuterVelX_DY" : ((2.0*y*y/mu2-apos)*exp(-1.0 *( x * x + y * y))),
           "SolutionOuterVelY_DX" : ((apos - 2*x*x / mu2)*exp(-1.0 *( x * x + y * y))),
           "SolutionOuterVelY_DY" : (-2.0*x*y/mu2*exp(-1.0 *( x * x + y * y))),
           "SolutionInnerVelX_DX" : (aneg*2*x*y*exp(-1.0 *( x * x + y * y))),
           "SolutionInnerVelX_DY" : (aneg*(2*y*y-1)*exp(-1.0 *( x * x + y * y))),
           "SolutionInnerVelY_DX" : (aneg*(1-2*x*x)*exp(-1.0 *( x * x + y * y))),
           "SolutionInnerVelY_DY" : (aneg*(-2)*x*y*exp(-1.0 *( x * x + y * y))),
           "SolutionInnerPressure" : (x*x*x),# + q),
           "SolutionOuterPressure" : (x*x*x - gammaf), #(pi*R*R/4.0*gammaf)),
}

coef_g = [ CoefficientFunction((functions["SourceInnerX"],functions["SourceInnerY"])),
           CoefficientFunction((functions["SourceOuterX"],functions["SourceOuterY"])) ]

vel_sol_neg = CoefficientFunction((functions["SolutionInnerVelX"],functions["SolutionInnerVelY"]))
vel_sol_pos = CoefficientFunction((functions["SolutionOuterVelX"],functions["SolutionOuterVelY"]))

grad_vel_sol_neg = CoefficientFunction((functions["SolutionInnerVelX_DX"],functions["SolutionInnerVelX_DY"],functions["SolutionInnerVelY_DX"],functions["SolutionInnerVelY_DY"]))
grad_vel_sol_pos = CoefficientFunction((functions["SolutionOuterVelX_DX"],functions["SolutionOuterVelX_DY"],functions["SolutionOuterVelY_DX"],functions["SolutionOuterVelY_DY"]))

pres_sol_neg = functions["SolutionInnerPressure"]
pres_sol_pos = functions["SolutionOuterPressure"]

# stabilization parameter for ghost-penalty
gamma_stab = 0.05
# stabilization parameter for Nitsche
lambda_nitsche  = 0.5 * (mu[0]+mu[1]) * 20 * order * order

levelset = functions["Levelset"]

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=10.5, discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lsetp1 = lsetmeshadap.lset_p1
ci = CutInfo(mesh,lsetp1)

Vh = H1(mesh, order=order, dirichlet=[1,2,3,4])
Vhx = XFESpace(Vh, ci)
VhG = FESpace([Vh, Vhx])
Qh = H1(mesh, order=order-1)
Qhx = XFESpace(Qh, ci)
QhG = FESpace([Qh, Qhx])
N = FESpace("number",mesh)

WhG = FESpace([VhG,VhG,QhG,N],flags = {"dgjumps" : True})

gfup = GridFunction(WhG)
gfu1 = gfup.components[0]
gfu2 = gfup.components[1]
gfp = gfup.components[2]
u1_coef = gfu1.components[0] + IfPos(lsetp1, pos(gfu1.components[1]), neg(gfu1.components[1]))
u2_coef = gfu2.components[0] + IfPos(lsetp1, pos(gfu2.components[1]), neg(gfu2.components[1]))

vel_neg = CoefficientFunction((gfu1.components[0]+neg(gfu1.components[1]),
                               gfu2.components[0]+neg(gfu2.components[1])))
vel_pos = CoefficientFunction((gfu1.components[0]+pos(gfu1.components[1]),
                               gfu2.components[0]+pos(gfu2.components[1])))
vel_coef = CoefficientFunction((u1_coef,u2_coef))

grad_vel_neg = CoefficientFunction((grad(gfu1.components[0])+neg_grad(gfu1.components[1]),
                                    grad(gfu2.components[0])+neg_grad(gfu2.components[1])))
grad_vel_pos = CoefficientFunction((grad(gfu1.components[0])+pos_grad(gfu1.components[1]),
                                    grad(gfu2.components[0])+pos_grad(gfu2.components[1])))


pres_neg = gfp.components[0]+neg(gfp.components[1])
pres_pos = gfp.components[0]+pos(gfp.components[1])
pres_coef = gfp.components[0] + IfPos(lsetp1, pos(gfp.components[1]), neg(gfp.components[1]))

gfp = gfup.components[2]

n_outer = specialcf.normal(mesh.dim)
h = specialcf.mesh_size   

lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
lset_pos = { "levelset" : lsetp1, "domain_type" : POS, "subdivlvl" : 0}
lset_if  = { "levelset" : lsetp1, "domain_type" : IF , "subdivlvl" : 0}

# element, facet and dof marking w.r.t. boundary approximation with lsetp1:
hasneg = ci.GetElementsOfType(HASNEG)
haspos = ci.GetElementsOfType(HASPOS)
hasif = ci.GetElementsOfType(IF)
ba_facets = [GetFacetsWithNeighborTypes(mesh,a=hasneg,b=hasif),GetFacetsWithNeighborTypes(mesh,a=haspos,b=hasif)]

n_lset = 1.0/Norm(grad(lsetp1)) * grad(lsetp1)
kappa = [CutRatioGF(ci),1.0-CutRatioGF(ci)]
           
a = BilinearForm(WhG,symmetric=False)
f = LinearForm(WhG)
            
(u1s,u1x),(u2s,u2x),(ps,px),n = WhG.TrialFunction()
(v1s,v1x),(v2s,v2x),(qs,qx),m = WhG.TestFunction()

u1 = [u1s + op(u1x) for op in [neg,pos]] 
u2 = [u2s + op(u2x) for op in [neg,pos]]
u = [CoefficientFunction((u1[i],u2[i])) for i in range(2)]
p = [ps + op(px) for op in [neg,pos]] 
dpdn = [grad(ps)*n_outer + op(px)*n_outer for op in [neg_grad,pos_grad]] 
dpdn_Other = [grad(ps.Other())*n_outer + op(px.Other())*n_outer for op in [neg_grad,pos_grad]] 
dpdn_jump = [dpdn[i] - dpdn_Other[i] for i in range(2)] 
gradu1 = [grad(u1s) + op(u1x) for op in [neg_grad,pos_grad]]
gradu2 = [grad(u2s) + op(u2x) for op in [neg_grad,pos_grad]]
divu = [gradu1[i][0]+gradu2[i][1] for i in range(2)]
Du = [CoefficientFunction((2*gradu1[i][0],gradu2[i][0]+gradu1[i][1],gradu2[i][0]+gradu1[i][1],2*gradu2[i][1]),dims = (2,2)) for i in range(2)]
sigmaupn = [ - mu[i] * (Du[i] * n_lset) + p[i] * n_lset for i in range(2)] 
average_flux_u = kappa[0] * sigmaupn[0] + kappa[1] * sigmaupn[1]

v1 = [v1s + op(v1x) for op in [neg,pos]] 
v2 = [v2s + op(v2x) for op in [neg,pos]] 
v = [CoefficientFunction((v1[i],v2[i])) for i in range(2)]
q = [qs + op(qx) for op in [neg,pos]] 
dqdn = [grad(qs)*n_outer + op(qx)*n_outer for op in [neg_grad,pos_grad]] 
dqdn_Other = [grad(qs.Other())*n_outer + op(qx.Other())*n_outer for op in [neg_grad,pos_grad]] 
dqdn_jump = [dqdn[i] - dqdn_Other[i] for i in range(2)] 
gradv1 = [grad(v1s) + op(v1x) for op in [neg_grad,pos_grad]]
gradv2 = [grad(v2s) + op(v2x) for op in [neg_grad,pos_grad]]
divv = [gradv1[i][0]+gradv2[i][1] for i in range(2)]
Dv = [CoefficientFunction((2*gradv1[i][0],gradv2[i][0]+gradv1[i][1],gradv2[i][0]+gradv1[i][1],2*gradv2[i][1]),dims = (2,2)) for i in range(2)]
sigmavqn = [ - mu[i] * Dv[i] * n_lset + q[i] * n_lset for i in range(2)]
average_flux_v = kappa[0] * sigmavqn[0] + kappa[1] * sigmavqn[1]
average_inv_v = - kappa[1] * v[0] - kappa[0] * v[1]

# Viscosity term
a += SymbolicBFI(lset_neg, form = 0.5 * mu1 * InnerProduct(Du[0],Dv[0]) )
a += SymbolicBFI(lset_pos, form = 0.5 * mu2 * InnerProduct(Du[1],Dv[1]) )
a += SymbolicBFI(lset_if,  form = InnerProduct( average_flux_u, v[0]-v[1]))
a += SymbolicBFI(lset_if,  form = InnerProduct( average_flux_v, u[0]-u[1]))
a += SymbolicBFI(lset_if,  form = lambda_nitsche/h * InnerProduct(u[0]-u[1],v[0]-v[1]))
# pressure stuff
a += SymbolicBFI(lset_neg, form = - divu[0] * q[0] - divv[0] * p[0]
                                  + n * q[0] + m * p[0])
a += SymbolicBFI(lset_pos, form = - divu[1] * q[1] - divv[1] * p[1])
#                                  + n * q[1] + m * p[1])

f += SymbolicLFI(lset_if,  form = gammaf * InnerProduct(average_inv_v,n_lset))
f += SymbolicLFI(lset_neg, form = coef_g[0] * v[0])
f += SymbolicLFI(lset_pos, form = coef_g[1] * v[1])


# # ghost penalty terms:
for i in range(2):
    a += SymbolicBFI(form = - gamma_stab * h*h*h* dpdn_jump[i]*dqdn_jump[i],VOL_or_BND = VOL, skeleton=True, definedonelements=ba_facets[i])

# apply mesh adaptation    
mesh.SetDeformation(deformation)

a.Assemble()
f.Assemble()

gfu1.components[0].Set(functions["SolutionOuterVelX"])
gfu2.components[0].Set(functions["SolutionOuterVelY"])
# Solve linear system
f.vec.data -= a.mat * gfup.vec
gfup.vec.data += a.mat.Inverse(WhG.FreeDofs()) * f.vec
              
# #measure the error

vl2error = sqrt( Integrate(lset_neg,InnerProduct(vel_neg-vel_sol_neg,vel_neg-vel_sol_neg),mesh=mesh)
                 +Integrate(lset_pos,InnerProduct(vel_pos-vel_sol_pos,vel_pos-vel_sol_pos),mesh=mesh) )
print("L2 Error of velocity: {0}".format(vl2error))

vh1error = sqrt( Integrate(lset_neg,InnerProduct(grad_vel_neg-grad_vel_sol_neg,grad_vel_neg-grad_vel_sol_neg),mesh=mesh)
                 +Integrate(lset_pos,InnerProduct(grad_vel_pos-grad_vel_sol_pos,grad_vel_pos-grad_vel_sol_pos),mesh=mesh) )
print("H1 Error of velocity: {0}".format(vh1error))

pl2error = sqrt(  Integrate(lset_neg,(pres_neg-pres_sol_neg)*(pres_neg-pres_sol_neg),mesh=mesh)
                 +Integrate(lset_pos,(pres_pos-pres_sol_pos)*(pres_pos-pres_sol_pos),mesh=mesh) )
print("L2 Error of pressure: {0}".format(pl2error))

# # unset mesh adaptation
mesh.UnsetDeformation()

# #visualization:
Draw(deformation,mesh,"deformation")
Draw(levelset,mesh,"levelset")
Draw(lsetp1,mesh,"lsetp1")
Draw(pres_coef,mesh,"pressure")
Draw(vel_coef,mesh,"vel")
