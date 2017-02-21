from math import pi
# ngsolve stuff
from ngsolve import *
# visualization stuff
from ngsolve.internal import *
# basic xfem functionality
from xfem import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

from netgen.geom2d import SplineGeometry

square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.2, quad_dominated=False))

r44 = (x*x*x*x+y*y*y*y)
r41 = sqrt(sqrt(x*x*x*x+y*y*y*y))
r4m3 = (1.0/(r41*r41*r41))
r66 = (x*x*x*x*x*x+y*y*y*y*y*y)
r63 = sqrt(r66)
r22 = (x*x+y*y)
r21 = sqrt(r22)

levelset = (sqrt(sqrt(x*x*x*x+y*y*y*y)) - 1.0)

solution = [1.0+pi/2.0-sqrt(2.0)*cos(pi/4.0*r44),pi/2.0*r41]
alpha = [1.0,2.0]
coef_f = [ (-1.0*sqrt(2.0)*pi*(pi*cos(pi/4*(r44))*(r66)+3*sin(pi/4*(r44))*(r22))),
          (-2.0*pi*3/2*(r4m3)*(-(r66)/(r44)+(r22))) ]

order = 2

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lsetp1 = lsetmeshadap.lset_p1

# extended FESpace 

Vh = H1(mesh, order=order, dirichlet=[1,2,3,4])
Vhx = XFESpace(Vh, mesh, lsetp1)
VhG = FESpace([Vh,Vhx])

# coefficients / parameters: 

n = 1.0/grad(lsetp1).Norm() * grad(lsetp1)
h = specialcf.mesh_size

kappa = kappa(mesh,lsetp1)

stab = 10*(alpha[1]+alpha[0])*(order+1)*order/h

# expressions of test and trial functions:

u_std, u_x = VhG.TrialFunction()
v_std, v_x = VhG.TestFunction()

u = [u_std + op(u_x) for op in [neg,pos]]
v = [v_std + op(v_x) for op in [neg,pos]]

gradu = [grad(u_std) + op(u_x) for op in [neg_grad,pos_grad]]
gradv = [grad(v_std) + op(v_x) for op in [neg_grad,pos_grad]]

average_flux_u = sum([- kappa[i] * alpha[i] * gradu[i] * n for i in [0,1]])
average_flux_v = sum([- kappa[i] * alpha[i] * gradv[i] * n for i in [0,1]])

# integration domains (and integration parameter "subdivlvl" and "force_intorder")

lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
lset_pos = { "levelset" : lsetp1, "domain_type" : POS, "subdivlvl" : 0}
lset_if  = { "levelset" : lsetp1, "domain_type" : IF , "subdivlvl" : 0}

# bilinear forms:

a = BilinearForm(VhG, symmetric = True, flags = { })
a += SymbolicBFI(levelset_domain = lset_neg, form = alpha[0] * gradu[0] * gradv[0])
a += SymbolicBFI(levelset_domain = lset_pos, form = alpha[1] * gradu[1] * gradv[1])
a += SymbolicBFI(levelset_domain = lset_if , form =     average_flux_u * (v[0]-v[1]))
a += SymbolicBFI(levelset_domain = lset_if , form =     average_flux_v * (u[0]-u[1]))
a += SymbolicBFI(levelset_domain = lset_if , form = stab * (u[0]-u[1]) * (v[0]-v[1]))

f = LinearForm(VhG)
f += SymbolicLFI(levelset_domain = lset_neg, form = coef_f[0] * v[0])
f += SymbolicLFI(levelset_domain = lset_pos, form = coef_f[1] * v[1])

gfu = GridFunction(VhG)

gfu.components[0].Set(solution[1], BND)

mesh.SetDeformation(deformation)

a.Assemble();
f.Assemble();

rhs = gfu.vec.CreateVector()
rhs.data = f.vec - a.mat * gfu.vec
update = gfu.vec.CreateVector()
update.data = a.mat.Inverse(VhG.FreeDofs()) * rhs;
gfu.vec.data += update

sol_coef = IfPos(lsetp1,solution[1],solution[0])
u_coef = gfu.components[0] + IfPos(lsetp1, pos(gfu.components[1]), neg(gfu.components[1]))
u = [gfu.components[0] + op(gfu.components[1]) for op in [neg,pos]]



err_sqr_coefs = [ (u[i] - solution[i])*(u[i] - solution[i]) for i in [0,1] ]

l2error = sqrt( Integrate( levelset_domain=lset_neg, cf=err_sqr_coefs[0], mesh=mesh,
                           order=2*order, heapsize=1000000)
                + Integrate(levelset_domain=lset_pos, cf=err_sqr_coefs[1], mesh=mesh,
                            order=2*order, heapsize=1000000))

print("L2 error : ",l2error)

mesh.UnsetDeformation()

Draw(lsetp1,mesh,"lsetp1")
Draw(lsetmeshadap.deform,mesh,"deformation")
Draw(u_coef,mesh,"u")
Draw(CoefficientFunction((lsetmeshadap.deform[0],lsetmeshadap.deform[1],u_coef)),mesh,"graph_of_u",sd=3)

visoptions.deformation = 1


