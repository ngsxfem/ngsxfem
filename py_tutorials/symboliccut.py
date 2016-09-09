from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *

from netgen.geom2d import SplineGeometry

# geometry

square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.7, quad_dominated=False))

levelset = sqrt(x*x+y*y) - 0.7
order = 1

lset_approx = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lset_approx)

# lset_approx = GridFunction(H1(mesh,order=order))
# lset_approx.Set(levelset)

# extended FESpace 

VhG = XStdFESpace(mesh, lset_approx, order=order, basetype="h1ho", dirichlet=[1,2,3,4])
gfu = GridFunction(VhG)

# coefficients / parameters: 

n = 1.0/sqrt(InnerProduct(grad(lset_approx),grad(lset_approx))) * grad(lset_approx)
h = specialcf.mesh_size

alpha_neg = 2.0
alpha_pos = 3.0
beta_neg = 1.0
beta_pos = 2.0

kappa_neg, kappa_pos = kappa(mesh,lset_approx)

stab = 10*(alpha_pos+alpha_neg)*order*order/h

# expressions of test and trial functions:

u_std, u_x = VhG.TrialFunction()
v_std, v_x = VhG.TestFunction()

u_pos = u_std + pos(u_x)
u_neg = u_std + neg(u_x)

v_pos = v_std + pos(v_x)
v_neg = v_std + neg(v_x)

gradu_pos = grad(u_std) + pos_grad(u_x)
gradu_neg = grad(u_std) + neg_grad(u_x)

gradv_pos = grad(v_std) + pos_grad(v_x)
gradv_neg = grad(v_std) + neg_grad(v_x)

betajump_u = beta_pos * u_pos - beta_neg * u_neg
betajump_v = beta_pos * v_pos - beta_neg * v_neg

average_flux_u = - kappa_pos * alpha_pos * gradu_pos * n - kappa_neg * alpha_neg * gradu_neg * n
average_flux_v = - kappa_pos * alpha_pos * gradv_pos * n - kappa_neg * alpha_neg * gradv_neg * n

# integration domains (and integration parameter "subdivlvl" and "force_intorder")

lset_neg = { "levelset" : lset_approx, "domain_type" : NEG, "subdivlvl" : 0}
lset_pos = { "levelset" : lset_approx, "domain_type" : POS, "subdivlvl" : 0}
lset_if  = { "levelset" : lset_approx, "domain_type" : IF , "subdivlvl" : 0}

# bilinear forms:

a = BilinearForm(VhG, symmetric = True, flags = { })
a += SymbolicBFI(levelset_domain = lset_neg, form = beta_neg * alpha_neg * gradu_neg * gradv_neg)
a += SymbolicBFI(levelset_domain = lset_pos, form = beta_pos * alpha_pos * gradu_pos * gradv_pos)
a += SymbolicBFI(levelset_domain = lset_if , form =  average_flux_u * betajump_v
                                                   + average_flux_v * betajump_u
                                                   + stab * betajump_u * betajump_v)
a.Assemble()


f = LinearForm(VhG)
f += SymbolicLFI(levelset_domain = lset_neg, form = 10 * v_neg)
f.Assemble();

gfu.components[0].Set(CoefficientFunction(0), boundary=True)

res = f.vec.CreateVector()
res.data = f.vec - a.mat * gfu.vec.data
gfu.vec.data += a.mat.Inverse(VhG.FreeDofs()) * res


u = gfu.components[0] + IfPos(lset_approx, pos(gfu.components[1]), neg(gfu.components[1]))

Draw(gfu.components[0],mesh,"u_std")
Draw(extend(gfu.components[1]),mesh,"u_x")
Draw(u,mesh,"u")



