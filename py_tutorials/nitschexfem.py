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

### Setting up discrete variational problem
VhG = XStdFESpace(mesh, lsetp1, order=order, basetype="h1ho", dirichlet=[1,2,3,4])

a = BilinearForm(VhG, symmetric = True, flags = { })
a += TwoDomainLaplaceIntegrator(alpha[0],alpha[1])

f = LinearForm(VhG)
f += TwoDomainSourceIntegrator(coef_f[0],coef_f[1])

nitsche_lhs, nitsche_jumpint, nitsche_fluxjumpint = XNitscheIntegrators(alpha, stab_param = 20.0)

a += nitsche_lhs

if nitsche_jumpint != None:
    f += nitsche_jumpint

if nitsche_fluxjumpint != None:
    f += nitsche_fluxjumpint

c = Preconditioner(a, type="local", flags= { "test" : True })

u = GridFunction(VhG)

u.components[0].Set(solution[1], BND)

mesh.SetDeformation(deformation)

a.Assemble();
f.Assemble();
c.Update();

solvea = CGSolver( mat=a.mat, pre=c.mat, complex=False, printrates=False, precision=1e-8, maxsteps=200000)
rhs = u.vec.CreateVector()
rhs.data = f.vec - a.mat * u.vec
update = u.vec.CreateVector()
update.data = solvea * rhs;
u.vec.data += update


#global last_num_its
#last_num_its = solvea.GetSteps()
mesh.UnsetDeformation()


sol_coef = IfPos(lsetp1,solution[1],solution[0])

# Draw(lsetp1,mesh,"lsetp1")
# # Draw(lsetmeshadap.deform,mesh,"deformation")
# Draw(u,mesh,"u")
# Draw(u-sol_coef,mesh,"err")

err_sqr_coefs = [ (u - solution[i])*(u - solution[i]) for i in [0,1] ]

l2error = sqrt(IntegrateOnWholeDomain(lsetp1, mesh, order=2*order,
                                      cf_neg = err_sqr_coefs[0], cf_pos = err_sqr_coefs[1]))
print("L2 error : ",l2error)
