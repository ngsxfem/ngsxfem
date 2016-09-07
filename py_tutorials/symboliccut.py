from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *

from netgen.geom2d import SplineGeometry

square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.8, quad_dominated=False))

levelset = sqrt(x*x+y*y) - 1.0
lset_approx = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lset_approx)

order = 1

VhG = XStdFESpace(mesh, lset_approx, order=order, basetype="h1ho", dirichlet=[])

u_std,u_x = VhG.TrialFunction()
v_std,v_x = VhG.TrialFunction()

gfu = GridFunction(VhG)

gfu.vec[:] = 1.0

gfu_pos = gfu.components[0] + pos(gfu.components[1])
gfu_neg = gfu.components[0] + neg(gfu.components[1])


u_pos = u_std + pos(u_x)
u_neg = u_std + neg(u_x)

v_pos = v_std + pos(v_x)
v_neg = v_std + neg(v_x)

a = BilinearForm(VhG, symmetric = True, flags = { })
# a += SymbolicCutBFI(lset=lset_approx,coef=u_neg*v_neg,domain_type=NEG)
# a += SymbolicCutBFI(lset=lset_approx,coef=u_pos*v_pos,domain_type=POS)
a += SymbolicBFI(u_std*v_std)
a += SymbolicBFI(extend(u_x)*extend(v_x))
a.Assemble()

f = LinearForm(VhG)
# f += TwoDomainSourceIntegrator(1,0)
f += SymbolicLFI(v_std)
f += SymbolicLFI(extend(v_x))
f.Assemble();

gfu.vec.data = a.mat.Inverse(VhG.FreeDofs()) * f.vec

# Draw(gfu.components[0]+pos(gfu.components[1]),mesh,"u")
Draw(gfu_pos,mesh,"u_pos")
Draw(gfu_neg,mesh,"u_neg")
Draw(gfu.components[0],mesh,"u_std")
Draw(extend(gfu.components[1]),mesh,"u_x")



