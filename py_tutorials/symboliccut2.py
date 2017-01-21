from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *

from netgen.geom2d import SplineGeometry

# geometry

square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))
mesh.Refine()
mesh.Refine()

levelset = sqrt(x*x+y*y) - 0.7
order = 1

lset_approx = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lset_approx)

# lset_approx = GridFunction(H1(mesh,order=order))
# lset_approx.Set(levelset)

# extended FESpace 

bulkfes = H1(mesh, order=order, dirichlet=[1,2,3,4])
VhG = FESpace([bulkfes,bulkfes])
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

u_neg, u_pos = VhG.TrialFunction()
v_neg, v_pos = VhG.TestFunction()

gradu_pos = grad(u_pos)
gradu_neg = grad(u_neg)

gradv_pos = grad(v_pos)
gradv_neg = grad(v_neg)

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

gfu.components[0].Set(CoefficientFunction(0),BND)
gfu.components[1].Set(CoefficientFunction(0),BND)

res = f.vec.CreateVector()
res.data = f.vec - a.mat * gfu.vec.data
gfu.vec.data += a.mat.Inverse(VhG.FreeDofs(),  inverse="sparsecholesky") * res


u = IfPos(lset_approx, gfu.components[1], gfu.components[0])

Draw(gfu.components[0],mesh,"u_neg")
Draw(gfu.components[1],mesh,"u_pos")
Draw(u,mesh,"u")



