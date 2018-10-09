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
            - alpha_i (u_xx + u_yy) =   f in subdomain i, i=1,2,
    interface conditions:
                                [u] =    0 on interface (continuity across the interface     ),
                   [-alpha · du/dn] =    0 on interface (conservation of the (diffusive) flux),
                                 u  =  u_D on domain boundary.

    The r.h.s. term f and the Dirichlet data u_D is chosen according to a manufactured solution
    which allows us to measure errors after the computation of a discrete solution.
    The coefficients alpha are domain-wise constants which are different in the two subdomains.

    discretization:
    ---------------
    Finite element space:
    For each subdomain we consider the space of piecewise linears restricted to the corresponding
    subdomain. 

    Variational formulation:
    We use a Nitsche formulation which involves averages of the fluxes and jumps of the solution
    across the interface. For the average we use the Hansbo-choice [1] where the average is adjusted
    to the local cut configuration in order to ensure stability of the resulting formulation.

    implementational aspects:
    ---------------
    Finite element space:
    The finite element space is a simple composition of two scalar piecewise linear finite element
    spaces w.r.t. the background mesh. However, in the variational formulation only these degrees of
    freedom that have some part in the corresponding subdomain are used and appear in the
    equations. Correspondingly, degrees of freedom in the active part of the mesh and unused degrees
    of freedoms are marked differently to obtain a consistent treatment and a regular stiffness
    matrix. 
    
    Geometry approximation:
    To approximate the implicitly described geometry we use the piecewise (multi-) linear
    interpolation of the level set function. For this geometry approximation (arbitrary order)
    accurate numerical integration is provided.

    linear systems:
    ---------------
    A (sparse) direct solver is applied to solve the arising linear systems.

    extensions:
    -----------
    * Instead of using the CutFEM-characterization of the space, one could use an
    XFEM-characterization, cf. nxfem.py. 
    * To obtain higher order accuracy (also w.r.t. the geometry approximatin) isoparametric unfitted
    methods can be used, cf. nxfem_higher_order.py.

    literature:
    -----------
    [1]: A.Hansbo, P.Hansbo, An unfitted finite element method, based on Nitsche's method, for
    elliptic interface problems, Comp. Meth. Appl. Mech. Eng., 2002

"""

# the constant pi
from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
# basic geometry features (for the background mesh)
from netgen.geom2d import SplineGeometry



# We generate the background mesh of the domain and use a simplicial triangulation
# To obtain a mesh with quadrilaterals use 'quad_dominated=True'
square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.2, quad_dominated=False))

# manufactured solution and corresponding r.h.s. data CoefficientFunctions:
r44 = (x*x*x*x+y*y*y*y)
r41 = sqrt(sqrt(x*x*x*x+y*y*y*y))
r4m3 = (1.0/(r41*r41*r41))
r66 = (x*x*x*x*x*x+y*y*y*y*y*y)
r63 = sqrt(r66)
r22 = (x*x+y*y)
r21 = sqrt(r22)
solution = [1.0+pi/2.0-sqrt(2.0)*cos(pi/4.0*r44),pi/2.0*r41]
coef_f = [ (-1.0*sqrt(2.0)*pi*(pi*cos(pi/4*(r44))*(r66)+3*sin(pi/4*(r44))*(r22))),
          (-2.0*pi*3/2*(r4m3)*(-(r66)/(r44)+(r22))) ]

# diffusion cofficients for the subdomains (NEG/POS):
alpha = [1.0,2.0]

# level set function of the domain (phi = ||x||_4 - 1) and its interpolation:
levelset = (sqrt(sqrt(x*x*x*x+y*y*y*y)) - 1.0)
lsetp1 = GridFunction(H1(mesh,order=1))
InterpolateToP1(levelset,lsetp1)

# Background FESpaces (used as CutFESpaces lateron):
Vh = H1(mesh, order=1, dirichlet=[1,2,3,4])
VhG = FESpace([Vh,Vh])
print("unknowns in background FESpace (2 x standard unknowns): ", VhG.ndof)

# Gathering information on cut elements:
#  * domain of (volume/boundary) element:
#    * NEG= only negative level set values
#    * POS= only positive level set values
#    * IF= cut element (negative and positive) level set values
#  * cut ratio:
#    If element is cut this describes the ratio between the measure of part in the negative domain
#    and the measure of the full element.
ci = CutInfo(mesh, lsetp1)

# Overwrite freedofs (degrees of freedoms that should be solved for) of VhG to mark only dofs that
# are involved in the cut problem. Use cut information of ci here:
hasneg = ci.GetElementsOfType(HASNEG)  # <- "hasneg": has (also) negative level set values
haspos = ci.GetElementsOfType(HASPOS)  # <- "haspos": has (also) positive level set values
freedofs = VhG.FreeDofs()
freedofs &= CompoundBitArray([GetDofsOfElements(Vh,hasneg),GetDofsOfElements(Vh,haspos)])

# coefficients / parameters:
n = 1.0/grad(lsetp1).Norm() * grad(lsetp1)
h = specialcf.mesh_size
# the cut ratio extracted from the cutinfo-class
kappa = (CutRatioGF(ci),1.0-CutRatioGF(ci))
# Nitsche stabilization parameter:
stab = 20*(alpha[1]+alpha[0])/h

# expressions of test and trial functions (u and v are tuples):
u = VhG.TrialFunction()
v = VhG.TestFunction()

gradu = [grad(ui) for ui in u]
gradv = [grad(vi) for vi in v]

average_flux_u = sum([- kappa[i] * alpha[i] * gradu[i] * n for i in [0,1]])
average_flux_v = sum([- kappa[i] * alpha[i] * gradv[i] * n for i in [0,1]])

# Integration domains for integration on negative/positive subdomains and on the interface:
# Here, the integration is (geometrically) exact if the "levelset"-argument is a piecewise
# (multi-)linear function. The integration order is chosen according to the arguments in the
# multilinear forms (but can be overwritten with "force_intorder" in the integration domain). If the
# "levelset"-argument is not a (multi-)linear function, you can use the "subdivlvl" argument to add
# additional refinement levels for the geometry approximation. 
lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
lset_pos = { "levelset" : lsetp1, "domain_type" : POS, "subdivlvl" : 0}
lset_if  = { "levelset" : lsetp1, "domain_type" : IF , "subdivlvl" : 0}

# bilinear forms:
a = BilinearForm(VhG, symmetric = True)
# l.h.s. domain integrals:
a += SymbolicBFI(levelset_domain = lset_neg, form = alpha[0] * gradu[0] * gradv[0])
a += SymbolicBFI(levelset_domain = lset_pos, form = alpha[1] * gradu[1] * gradv[1])
# Nitsche integrals:
a += SymbolicBFI(levelset_domain = lset_if , form =       average_flux_u * (v[0]-v[1])
                                                    +     average_flux_v * (u[0]-u[1])
                                                    + stab * (u[0]-u[1]) * (v[0]-v[1]))

f = LinearForm(VhG)
# r.h.s. domain integrals:
f += SymbolicLFI(levelset_domain = lset_neg, form = coef_f[0] * v[0])
f += SymbolicLFI(levelset_domain = lset_pos, form = coef_f[1] * v[1])

# solution vector
gfu = GridFunction(VhG)

# setting domain boundary conditions:
gfu.components[1].Set(solution[1], BND)

# setting up matrix and vector
a.Assemble()
f.Assemble()

# homogenization of boundary data and solution of linear system
rhs = gfu.vec.CreateVector()
rhs.data = f.vec - a.mat * gfu.vec
update = gfu.vec.CreateVector()
update.data = a.mat.Inverse(freedofs) * rhs;
gfu.vec.data += update

# visualization of (discrete) solution: Wherever (interpolated) level set function is negative
# visualize the first component, where it is positive visualize the second component
u_coef = IfPos(lsetp1, gfu.components[1], gfu.components[0])

# visualize levelset, interpolated levelset and discrete solution:
# (Note that the visualization does not respect the discontinuities. They are smeared out. To see
#  kinks or jumps more clearly increase the subdivision option in the visualization.)
import sys
if not hasattr(sys, 'argv') or len(sys.argv) == 1 or sys.argv[1] != "testmode":
    Draw(levelset,mesh,"levelset")
    Draw(lsetp1,mesh,"levelset_P1")
    Draw(u_coef,mesh,"u")

# Error coefficients:
err_sqr_coefs = [(gfu.components[i]-solution[i])*(gfu.components[i]-solution[i]) for i in [0,1] ]

# Computation of L2 error:
l2error = sqrt(   Integrate( levelset_domain=lset_neg, cf=err_sqr_coefs[0], mesh=mesh, order=2)
                + Integrate( levelset_domain=lset_pos, cf=err_sqr_coefs[1], mesh=mesh, order=2))

print("L2 error : ",l2error)
