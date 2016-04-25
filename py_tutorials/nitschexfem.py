from math import pi
from time import sleep
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *


from netgen.geom2d import SplineGeometry
square = SplineGeometry()
square.AddRectangle([-1,-1],[1,1],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.1, quad_dominated=False))

levelset = sqrt(x*x+y*y)-0.7

sleep(1)

solution = sin(pi*z)
alpha = [1.0,2.0]
coef_f = [x*x,y*y]
order = 2

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)

### Setting up discrete variational problem

VhG = XStdFESpace(mesh, levelset, order=1, basetype="h1ho", dirichlet=[1,2,3,4])

a = BilinearForm(VhG, symmetric = True, flags = { })
a += TwoDomainLaplaceIntegrator(alpha[0],alpha[1])

f = LinearForm(VhG)
f += TwoDomainSourceIntegrator(coef_f[0],coef_f[1])

nitsche_lhs, nitsche_jumpint, nitsche_fluxjumpint = XNitscheIntegrators(alpha,jump=1.0,fluxjump=1.0)

a += nitsche_lhs

if nitsche_jumpint != None:
    f += nitsche_jumpint

if nitsche_fluxjumpint != None:
    f += nitsche_fluxjumpint

c = Preconditioner(a, type="local", flags= { "test" : True })

u = GridFunction(VhG)


deformation = lsetmeshadap.CalcDeformation(levelset)
mesh.SetDeformation(deformation)

# # VhG.Update()
# # u.Update()

a.Assemble();
f.Assemble();
c.Update();

solvea = CGSolver( mat=a.mat, pre=c.mat, complex=False, printrates=False, precision=1e-8, maxsteps=200000)
u.vec.data = solvea * f.vec;

global last_num_its
last_num_its = solvea.GetSteps()
mesh.UnsetDeformation()

Draw(lsetmeshadap.lset_p1,mesh,"lsetp1")
Draw(lsetmeshadap.deform,mesh,"deformation")
Draw(u,mesh,"u")

