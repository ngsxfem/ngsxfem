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
for i in range(6):
    mesh.Refine()


mu1 = 1.0
mu2 = 10.0
mu = [mu1,mu2]

R = 2.0/3.0
aneg = 1.0/mu1
apos = 1.0/mu2 + (1.0/mu1 - 1.0/mu2)*exp(x*x+y*y-R*R)
gammaf = 0.5
q = gammaf - pi*R*R/4.0*gammaf

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
           "SolutionInnerPressure" : (x*x*x + q),
           "SolutionOuterPressure" : (x*x*x - (pi*R*R/4.0*gammaf)),
}

coef_g = [ CoefficientFunction((functions["SourceInnerX"],functions["SourceInnerY"])),
           CoefficientFunction((functions["SourceOuterX"],functions["SourceOuterY"])) ]

vel_sol_neg = CoefficientFunction((functions["SolutionInnerVelX"],functions["SolutionInnerVelY"]))
vel_sol_pos = CoefficientFunction((functions["SolutionOuterVelX"],functions["SolutionOuterVelY"]))

order = 2

# stabilization parameter for ghost-penalty
gamma_stab = 0.05
# stabilization parameter for Nitsche
lambda_nitsche  = 0.5 * (mu[0]+mu[1]) * 20 * order * order

levelset = functions["Levelset"]

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.5, discontinuous_qn=True)
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

WhG = FESpace([VhG,VhG,QhG,N],flags = {"dgjumps":True})

gfup = GridFunction(WhG)
gfu1 = gfup.components[0]
gfu2 = gfup.components[1]
u1_coef = gfu1.components[0] + IfPos(lsetp1, pos(gfu1.components[1]), neg(gfu1.components[1]))
u2_coef = gfu2.components[0] + IfPos(lsetp1, pos(gfu2.components[1]), neg(gfu2.components[1]))

vel_neg = CoefficientFunction((gfu1.components[0]+neg(gfu1.components[1]),
                               gfu2.components[0]+neg(gfu2.components[1])))
vel_pos = CoefficientFunction((gfu1.components[0]+pos(gfu1.components[1]),
                               gfu2.components[0]+pos(gfu2.components[1])))
vel_coef = CoefficientFunction((u1_coef,u2_coef))

gfp = gfup.components[2]

n_outer = specialcf.normal(mesh.dim)
h = specialcf.mesh_size   

lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
lset_pos = { "levelset" : lsetp1, "domain_type" : POS, "subdivlvl" : 0}
lset_if  = { "levelset" : lsetp1, "domain_type" : IF , "subdivlvl" : 0}

# element, facet and dof marking w.r.t. boundary approximation with lsetp1:
hasneg = BitArray(ci.GetElementsOfType(NEG))
hasneg |= ci.GetElementsOfType(IF)
haspos = BitArray(ci.GetElementsOfType(POS))
haspos |= ci.GetElementsOfType(IF)
hasif = BitArray(ci.GetElementsOfType(IF))
hasif |= ci.GetElementsOfType(IF) 
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
a += SymbolicBFI(lset_if,  form = lambda_nitsche * InnerProduct(u[0]-u[1],v[0]-v[1]))
# pressure stuff
a += SymbolicBFI(lset_neg, form = - divu[0] * q[0] - divv[0] * p[0]
                                  + n * q[0] + m * p[0])
a += SymbolicBFI(lset_pos, form = - divu[1] * q[1] - divv[1] * p[1]
                                  + n * q[1] + m * p[1])

f += SymbolicLFI(lset_if,  form = gammaf * InnerProduct(average_inv_v,n_lset))
f += SymbolicLFI(lset_neg, form = coef_g[0] * v[0])
f += SymbolicLFI(lset_pos, form = coef_g[1] * v[1])


# # ghost penalty terms:
for i in range(2):
    bfi_gp = SymbolicBFI(form = - gamma_stab * h*h*h* dpdn_jump[i]*dqdn_jump[i],VOL_or_BND = VOL, skeleton=True)
    bfi_gp.SetDefinedOnElements(ba_facets[i])
    a += bfi_gp

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
l2error = sqrt( Integrate(lset_neg,InnerProduct(vel_neg-vel_sol_neg,vel_neg-vel_sol_neg),mesh=mesh)
               +Integrate(lset_pos,InnerProduct(vel_pos-vel_sol_pos,vel_pos-vel_sol_pos),mesh=mesh) )
print("L2 Error: {0}".format(l2error))

# # unset mesh adaptation
mesh.UnsetDeformation()

# #visualization:

Draw(deformation,mesh,"deformation")
Draw(levelset,mesh,"levelset")
Draw(lsetp1,mesh,"lsetp1")
Draw(vel_coef,mesh,"vel")

