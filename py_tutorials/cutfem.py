"""
In this example we solve a scalar *unfitted* interface problem. As a discretization method we use a
level set based geometry description and a Cut Finite element method with a Nitsche formulation to
impose the interface conditions. 

    domain: 
    -------

    The domain is [-1.5,1.5]^2 while the interface is described by a level set function ( unit ball
    in the ||·||_4 norm ). In the discretization the level set function is approximated with a
    piecewise linear interpolation. 

    PDE problem:
    ------------
    domain equations:
      - alpha_i u_xx - alpha_i u_yy = f in subdomain i, i=1,2
    interface conditions:
                                [u] = 0 on interface (continuity across the interface     )
                   [-alpha · du/dn] = 0 on interface (conservation of the (diffusive) flux)

    discretization:
    ---------------
    Finite element space:
    For each subdomain we consider the space of piecewise linears restricted to the corresponding
    subdomain. 

    Variational formulation:
    We use a Nitsche formulation which involves averages of the fluxes and jumps of the solution
    across the interface. For the average we use the Hansbo-choice where the average is adjusted to
    the local cut configuration in order to ensure stability of the resulting formulation.
    
    ...to be continued...

    implementation:
    ---------------
    Finite element space:
    The finite element space is a simple composition of two scalar piecewise linear finite element
    spaces. However, in the variational formulation only these degrees of freedom that have some
    part in the corresponding subdomain are used and appear in the equations. Many unused dofs
    appear and have to be marked to finally obtain a regular matrix.

    ...to be continued...


"""

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
mesh = Mesh (square.GenerateMesh(maxh=0.2, quad_dominated=False))

r44 = (x*x*x*x+y*y*y*y)
r41 = sqrt(sqrt(x*x*x*x+y*y*y*y))
r4m3 = (1.0/(r41*r41*r41))
r66 = (x*x*x*x*x*x+y*y*y*y*y*y)
r63 = sqrt(r66)
r22 = (x*x+y*y)
r21 = sqrt(r22)

levelset = (sqrt(sqrt(x*x*x*x+y*y*y*y)) - 1.0)
lsetp1 = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lsetp1)

solution = [1.0+pi/2.0-sqrt(2.0)*cos(pi/4.0*r44),pi/2.0*r41]
alpha = [1.0,2.0]
coef_f = [ (-1.0*sqrt(2.0)*pi*(pi*cos(pi/4*(r44))*(r66)+3*sin(pi/4*(r44))*(r22))),
          (-2.0*pi*3/2*(r4m3)*(-(r66)/(r44)+(r22))) ]

order = 1

# Cut FESpaces:

Vh = H1(mesh, order=order, dirichlet=[1,2,3,4])
VhG = FESpace([Vh,Vh])
print("unknowns in CutFESpace (2 x standard unknowns): ", VhG.ndof)

# overwrite freedofs of VhG to mark only dofs that are involved in the cut problem
ci = CutInfo(mesh, lsetp1)
hasneg = BitArray(ci.GetElementsOfType(NEG))
hasneg |= ci.GetElementsOfType(IF)
haspos = BitArray(ci.GetElementsOfType(POS))
haspos |= ci.GetElementsOfType(IF)
freedofs = VhG.FreeDofs()
freedofs &= CompoundBitArray([GetDofsOfElements(Vh,hasneg),GetDofsOfElements(Vh,haspos)])

# coefficients / parameters:
n = 1.0/grad(lsetp1).Norm() * grad(lsetp1)
h = specialcf.mesh_size

kappa = kappa(mesh,lsetp1)

stab = 10*(alpha[1]+alpha[0])*(order+1)*order/h

# expressions of test and trial functions:

u = VhG.TrialFunction()
v = VhG.TestFunction()

gradu = [grad(ui) for ui in u]
gradv = [grad(vi) for vi in v]

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

gfu.components[1].Set(solution[1], BND)

a.Assemble();
f.Assemble();
rhs = gfu.vec.CreateVector()
rhs.data = f.vec - a.mat * gfu.vec
update = gfu.vec.CreateVector()
update.data = a.mat.Inverse(freedofs) * rhs;
gfu.vec.data += update

sol_coef = IfPos(lsetp1,solution[1],solution[0])
u_coef = IfPos(lsetp1, gfu.components[1], gfu.components[0])
u = [gfu.components[i] for i in [0,1]]

Draw(levelset,mesh,"levelset")
Draw(lsetp1,mesh,"levelset_P1")
Draw(u_coef,mesh,"u")

err_sqr_coefs = [ (u[i] - solution[i])*(u[i] - solution[i]) for i in [0,1] ]

l2error = sqrt( Integrate( levelset_domain=lset_neg, cf=err_sqr_coefs[0], mesh=mesh,
                           order=2*order, heapsize=1000000)
                + Integrate(levelset_domain=lset_pos, cf=err_sqr_coefs[1], mesh=mesh,
                            order=2*order, heapsize=1000000))

print("L2 error : ",l2error)
