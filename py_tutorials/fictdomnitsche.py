from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *

from netgen.geom2d import SplineGeometry
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

# geometry

square = SplineGeometry()
square.AddRectangle([-1,-1],[1,1],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.125, quad_dominated=False))

sol=y*y

levelset = sqrt(x*x+y*y) - 0.7
order = 2

lset_approx = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lset_approx)

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)

# lset_approx = GridFunction(H1(mesh,order=order))
# lset_approx.Set(levelset)

# extended FESpace 

Vh = H1(mesh, order=order, dirichlet=[1,2,3,4], flags = {"dgjumps" : True})
gfu = GridFunction(Vh)

# coefficients / parameters: 

n = 1.0/sqrt(InnerProduct(grad(lset_approx),grad(lset_approx))) * grad(lset_approx)
h = specialcf.mesh_size

u = Vh.TrialFunction()
v = Vh.TestFunction()

print("Vh.ndof:",Vh.ndof)

# alpha_neg = 2.0
# alpha_pos = 3.0
# beta_neg = 1.0
# beta_pos = 2.0

# kappa_neg, kappa_pos = kappa(mesh,lset_approx)

# stab = 10*(alpha_pos+alpha_neg)*order*order/h

# # expressions of test and trial functions:

# u_neg, u_pos = VhG.TrialFunction()
# v_neg, v_pos = VhG.TestFunction()

# gradu_pos = grad(u_pos)
# gradu_neg = grad(u_neg)

# gradv_pos = grad(v_pos)
# gradv_neg = grad(v_neg)

# betajump_u = beta_pos * u_pos - beta_neg * u_neg
# betajump_v = beta_pos * v_pos - beta_neg * v_neg

# average_flux_u = - kappa_pos * alpha_pos * gradu_pos * n - kappa_neg * alpha_neg * gradu_neg * n
# average_flux_v = - kappa_pos * alpha_pos * gradv_pos * n - kappa_neg * alpha_neg * gradv_neg * n

# integration domains (and integration parameter "subdivlvl" and "force_intorder")

# lset_neg = { "levelset" : lset_approx, "domain_type" : NEG, "subdivlvl" : 0}
# lset_pos = { "levelset" : lset_approx, "domain_type" : POS, "subdivlvl" : 0}
# lset_if  = { "levelset" : lset_approx, "domain_type" : IF , "subdivlvl" : 0}

# bilinear forms:
# for el in Vh.Elements():
#     print(el.dofs)

a = BilinearForm(Vh, symmetric = False, flags = { })
# a += SymbolicBFI(levelset_domain = lset_neg, form = beta_neg * alpha_neg * gradu_neg * gradv_neg)
# a += SymbolicBFI(levelset_domain = lset_pos, form = beta_pos * alpha_pos * gradu_pos * gradv_pos)
# a += SymbolicBFI(levelset_domain = lset_if , form =  average_flux_u * betajump_v
#                                                    + average_flux_v * betajump_u
#                                                    + stab * betajump_u * betajump_v)

def dnjump(u,order):
    if order%2==0:
        return dn(u,order) - dn(u.Other(),order)
    else:
        return dn(u,order) + dn(u.Other(),order)

factors = [1.0/h, h, h*h*h, h*h*h*h*h, h*h*h*h*h*h*h]
# a += SymbolicBFI( u * v ) 

for i in range(1,order+1):
    a += SymbolicBFI( factors[i] * dnjump(u,i) * dnjump(v,i), skeleton=True )

deformation = lsetmeshadap.CalcDeformation(levelset)
# mesh.SetDeformation(deformation)

a.Assemble()

f = LinearForm(Vh)
# f += SymbolicLFI(sol*v)
f.Assemble();

gfu.Set(sol, boundary=True)

# gfu.components[0].Set(CoefficientFunction(0), boundary=True)
# gfu.components[1].Set(CoefficientFunction(0), boundary=True)

res = f.vec.CreateVector()
res.data = f.vec - a.mat * gfu.vec.data
gfu.vec.data += a.mat.Inverse(Vh.FreeDofs()) * res

# u = IfPos(lset_approx, gfu.components[1], gfu.components[0])

# Draw(gfu.components[0],mesh,"u_neg")
# Draw(gfu.components[1],mesh,"u_pos")
# gfu.Set(y*y)
# print(gfu.vec)

# gfu.vec[:] = 0.0
# gfu.vec[8] = 1.0

# Lh = L2(mesh, order=order, dirichlet=[1,2,3,4], flags = {"dgjumps" : True})
# gfdudn = GridFunction(Lh)
# gfdudn.Set(-grad(gfu)[0]/sqrt(2.0)+grad(gfu)[1]/sqrt(2.0))
# gfdudn2 = GridFunction(Lh)
# gfdudn2.Set(-grad(gfdudn)[0]/sqrt(2.0)+grad(gfdudn)[1]/sqrt(2.0))

Draw(deformation,mesh,"deformation")

Draw(gfu,mesh,"u")
# Draw(gfdudn,mesh,"gfdudn")
# Draw(gfdudn2,mesh,"gfdudn2")

# Draw(dn(gfu,1),mesh,"ddu")



print(sqrt(Integrate_old((sol-gfu)*(sol-gfu),mesh,order=2*order)))

mesh.UnsetDeformation()
