from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

from netgen.geom2d import SplineGeometry
square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.1, quad_dominated=False))

levelset = sqrt(x*x+y*y) - 1.0
# level set approximation that is used in the implementation:
lset_approx = levelset

# diffusion coefficients
alpha = [2.0,1.0]
# source terms w.r.t. domain 1/2
coef_f = [1.0,0.0]
# r.h.s. g for interface condition [u]=g
jump = 1.0
# r.h.s. h for interface condition [alpha du/dn]=h
fluxjump = None
# stabilization parameter of Nitsche discretization
stab_param = 20.0
# polynomial degree
order = 3
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lset_approx = lsetmeshadap.lset_p1

VhG = XStdFESpace(mesh, lset_approx, order=order, basetype="h1ho", dirichlet=[1,2,3,4])

a = BilinearForm(VhG, symmetric = True, flags = { })
a += TwoDomainLaplaceIntegrator(alpha[0],alpha[1])

f = LinearForm(VhG)
f += TwoDomainSourceIntegrator(coef_f[0],coef_f[1])

nitsche_lhs, nitsche_jumpint, nitsche_fluxjumpint = XNitscheIntegrators(alpha, stab_param = stab_param, 
                                                                        jump=jump, fluxjump=fluxjump)

a += nitsche_lhs

if nitsche_jumpint != None:
    f += nitsche_jumpint

if nitsche_fluxjumpint != None:
    f += nitsche_fluxjumpint

u = GridFunction(VhG)

u.components[0].Set(sin(pi*(x+y)), BND)

mesh.SetDeformation(deformation)

a.Assemble();
f.Assemble();

mesh.UnsetDeformation()

rhs = u.vec.CreateVector()
rhs.data = f.vec - a.mat * u.vec
u.vec.data += a.mat.Inverse(VhG.FreeDofs()) * rhs

#Draw(lset_approx,mesh,"lset_approx")
#Draw(lsetmeshadap.deform,mesh,"deformation")
#Draw(u,mesh,"u")
