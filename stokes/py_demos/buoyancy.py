from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *
# For Stokes-FESpace and Stokes-Integrators (convenience)
from xfem.stokes import *

from netgen.geom2d import SplineGeometry
square = SplineGeometry()
square.AddRectangle([-1,-1],[1,1],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.1, quad_dominated=False))

levelset = sqrt(x*x+y*y) - 1.0
# level set approximation that is used in the implementation:
lset_approx = levelset

from netgen.geom2d import SplineGeometry
square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.1, quad_dominated=False))

levelset = sqrt(x*x+y*y) - 1.0
# level set approximation that is used in the implementation:
lset_approx = levelset

# diffusion coefficients (two domains)
mu = [2.0,1.0]
# density coefficients (two domains)
rho = [10.0,1.0]
# gravity constant (vector)
g = [0,10.0]
# interface force (e.g. surface tension)
gammaf = 1.0
# stabilization parameter of Nitsche discretization
stab_param = 20.0
# polynomial degree
order = 3
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lset_approx = lsetmeshadap.lset_p1


Vh = XStokesFESpace(mesh, order=order, levelset=lsetmeshadap.lset_p1, dirichlet=[1,2,3,4])

a = BilinearForm(Vh, symmetric = True, flags = { })
a += TwoDomainStokesIntegrator(mu[0], mu[1])
nitsche_a, nitsche_f = NitscheStokesIntegrators(mu[0], mu[1],
                                                lamb=stab_param,
                                                gammaf=gammaf)
a += nitsche_a

f = LinearForm(Vh)
f += nitsche_f
f.components[0] += TwoDomainSourceIntegrator(rho[0] * g[0],
                                             rho[1] * g[0])
f.components[1] += TwoDomainSourceIntegrator(rho[0] * g[1],
                                             rho[1] * g[1])

# uvp.components[0].components[0].Set(dir_u, boundary=True, definedon=mesh.Boundaries("dirbound"))
# uvp.components[1].components[0].Set(dir_v, boundary=True, definedon=mesh.Boundaries("dirbound"))

uvp = GridFunction(Vh)

mesh.SetDeformation(deformation)

a.Assemble();
f.Assemble();

mesh.UnsetDeformation()

f.vec.data -= a.mat * uvp.vec
uvp.vec.data += a.mat.Inverse(Vh.FreeDofs()) * f.vec

Draw(levelset,mesh,"levelset")
Draw(lsetmeshadap.deform,mesh,"deformation")
Draw(lsetmeshadap.lset_p1,mesh,"levelset(P1)")
Draw(CoefficientFunction (uvp.components[0:2]),mesh,"velocity")
Draw(CoefficientFunction (uvp.components[2]),mesh,"pressure")
