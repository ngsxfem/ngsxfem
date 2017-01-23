from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *


from netgen.geom2d import SplineGeometry

order = 2
threshold = 0.04 # curvate/resolution driven refinements (for trigs)
quad = True

#def get_l2error(order, n_ref, deform):
square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=10, quad_dominated=quad))
mesh.Refine()
mesh.Refine()

r44 = (x*x*x*x+y*y*y*y)
r41 = sqrt(sqrt(x*x*x*x+y*y*y*y))
r4m3 = (1.0/(r41*r41*r41))
r66 = (x*x*x*x*x*x+y*y*y*y*y*y)
r63 = sqrt(r66)
r22 = (x*x+y*y)
r21 = sqrt(r22)

levelset = (sqrt(sqrt(x*x*x*x+y*y*y*y)) - 1.0)

if not quad:
    for i in range(8):
        lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=threshold, discontinuous_qn=True)
        lsetmeshadap.CalcDeformation(levelset)
        # not working for quads:
        lsetmeshadap.MarkForRefinement(levelset,refine_threshold=threshold,absolute=False)
        mesh.Refine()

for i in range(3):
    mesh.Refine()


solution = [1.0+pi/2.0-sqrt(2.0)*cos(pi/4.0*r44),pi/2.0*r41]
alpha = [1.0,2.0]
coef_f = [ (-1.0*sqrt(2.0)*pi*(pi*cos(pi/4*(r44))*(r66)+3*sin(pi/4*(r44))*(r22))),
          (-2.0*pi*3/2*(r4m3)*(-(r66)/(r44)+(r22))) ]

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=threshold, discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lsetp1 = lsetmeshadap.lset_p1

# extended FESpace 

bulkfes = H1(mesh, order=order, dirichlet=[1,2,3,4])
VhG = FESpace([bulkfes,bulkfes])
#VhG = XStdFESpace(mesh, lsetp1, order=order, basetype="h1ho", dirichlet=[1,2,3,4])
gfu = GridFunction(VhG)

# coefficients / parameters: 

n = 1.0/sqrt(InnerProduct(grad(lsetp1),grad(lsetp1))) * grad(lsetp1)
h = specialcf.mesh_size

kappa_neg, kappa_pos = kappa(mesh,lsetp1)

stab = 20*(alpha[1]+alpha[0])/2.0*(order+1)*order/h

# expressions of test and trial functions:

u_neg, u_pos = VhG.TrialFunction()
v_neg, v_pos = VhG.TestFunction()

gradu_pos = grad(u_pos)
gradu_neg = grad(u_neg)

gradv_pos = grad(v_pos)
gradv_neg = grad(v_neg)

jump_u = -u_pos + u_neg
jump_v = -v_pos + v_neg

average_flux_u = - kappa_pos * alpha[1] * gradu_pos * n - kappa_neg * alpha[0] * gradu_neg * n
average_flux_v = - kappa_pos * alpha[1] * gradv_pos * n - kappa_neg * alpha[0] * gradv_neg * n

# integration domains (and integration parameter "subdivlvl" and "force_intorder")

lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
lset_pos = { "levelset" : lsetp1, "domain_type" : POS, "subdivlvl" : 0}
lset_if  = { "levelset" : lsetp1, "domain_type" : IF , "subdivlvl" : 0}

# bilinear forms:

a = BilinearForm(VhG, symmetric = True, flags = { })
# a += SymbolicBFI(levelset_domain = lset_neg, form = u_neg * v_neg)
# a += SymbolicBFI(levelset_domain = lset_pos, form = u_pos * v_pos)
a += SymbolicBFI(levelset_domain = lset_neg, form = alpha[0] * gradu_neg * gradv_neg)
a += SymbolicBFI(levelset_domain = lset_pos, form = alpha[1] * gradu_pos * gradv_pos)
a += SymbolicBFI(levelset_domain = lset_if , form =  average_flux_u * jump_v
                                                + average_flux_v * jump_u
                                                + stab * jump_u * jump_v)


f = LinearForm(VhG)
f += SymbolicLFI(levelset_domain = lset_neg, form = coef_f[0] * v_neg)
f += SymbolicLFI(levelset_domain = lset_pos, form = coef_f[1] * v_pos)

#c = Preconditioner(a, type="local", flags= { "test" : True })
c = Preconditioner(a, type="direct", flags= { "inverse" : "pardiso" })

u = GridFunction(VhG)

u.components[1].Set(solution[1], BND)

#if deform:
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


sol_coef = IfPos(lsetp1,solution[1],solution[0])
u_coef = IfPos(lsetp1,u.components[1],u.components[0])

Draw(lsetp1,mesh,"lsetp1")
Draw(lsetmeshadap.deform,mesh,"deformation")
Draw(u_coef,mesh,"u")
Draw(CoefficientFunction((lsetmeshadap.deform[0],lsetmeshadap.deform[1],u_coef)),mesh,"evelation_u",sd=4)

err_sqr_coefs = [ (u.components[i] - solution[i])*(u.components[i] - solution[i]) for i in [0,1] ]

l2error = sqrt(NewIntegrateX(lset=lsetp1,mesh=mesh,cf=err_sqr_coefs[0],order=2*order,domain_type=NEG,heapsize=1000000, use_saye=False) + NewIntegrateX(lset=lsetp1,mesh=mesh,cf=err_sqr_coefs[1],order=2*order,domain_type=POS,heapsize=1000000, use_saye=False))

mesh.UnsetDeformation()
print("L2 error : ",l2error)

