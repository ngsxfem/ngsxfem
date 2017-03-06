from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *


from netgen.csg import OrthoBrick
from netgen.csg import CSGeometry, OrthoBrick, Pnt
cube = OrthoBrick( Pnt(-1.5,-1.5,-1.5), Pnt(1.5,1.5,1.5) ).bc(1)
geom = CSGeometry()
geom.Add (cube)
mesh = Mesh (geom.GenerateMesh(maxh=10.1, quad_dominated=False))
# mesh.Refine()
# mesh.Refine()
# mesh.Refine()
# mesh.Refine()

levelset = sqrt(x*x+y*y+z*z) - 1.0
d = CoefficientFunction(1)
c = 1
sol = sin(pi*z)
rhs = sin(pi*z)*(d*pi*pi*(1-z*z)+c)+d*cos(pi*z)*2*pi*z
# sol = 1
# rhs = 1
order = 1

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lset_approx = lsetmeshadap.lset_p1

# extended FESpace 

VhG = H1(mesh, order=order)
gfu = GridFunction(VhG)

# coefficients / parameters: 

n = 1.0/sqrt(InnerProduct(grad(lset_approx),grad(lset_approx))) * grad(lset_approx)
h = specialcf.mesh_size

u = VhG.TrialFunction()
v = VhG.TestFunction()

lset_if  = { "levelset" : lset_approx, "domain_type" : IF , "subdivlvl" : 0}

# bilinear forms:

a = BilinearForm(VhG, symmetric = True, flags = { })
a += SymbolicBFI(levelset_domain = lset_if, form = d * InnerProduct(grad(u),grad(v)) )
a += SymbolicBFI(levelset_domain = lset_if, form = c * u * v)
# a += SymbolicBFI(form = (grad(u)*n) * (grad(v)*n))

f = LinearForm(VhG)
f += SymbolicLFI(levelset_domain = lset_if, form = rhs * v)

c = Preconditioner(a, type="direct", flags= { "inverse" : "pardiso" })

u = GridFunction(VhG)


mesh.SetDeformation(deformation)

a.Assemble();
f.Assemble();
c.Update();

solvea = CGSolver( mat=a.mat, pre=c.mat, complex=False,
                   printrates=True, precision=1e-8, maxsteps=200000)
u.vec.data = solvea * f.vec


#global last_num_its
#last_num_its = solvea.GetSteps()


Draw(lset_approx,mesh,"lsetp1")
Draw(u,mesh,"u")
Draw(deformation,mesh,"deformation")

l2error = sqrt(NewIntegrateX(lset=lset_approx,mesh=mesh,cf=(u-sol)*(u-sol),order=order,domain_type=IF,heapsize=1000000, use_saye=False))

print("L2 error : ", l2error)

mesh.UnsetDeformation()
