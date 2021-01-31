"""
In this example we solve an *unfitted* Stokes interface problem with a
high order isoparametric unfitted discretisation method. This file is
based on the discretisations in cutfem.py, nxfem.py and
nxfem_higher_order.py.

Domain:
-------
The domain is [-1,1]^2 while the interface is described by a level set
function ( circle with radius R=2/3).

PDE problem:
------------
domain equations for velocity (u1, u2) and pressure p:
 mu_i div( sigma(u_i,p_i) ) = rho_i g  in sub-domain i, i=1,2,
                   div(u_i) = 0        in sub-domain i, i=1,2,
interface conditions:
                        [u] = 0        on interface (continuity across
                                                     the interface),
             [sigma(u,p)·n] = f        on interface (conservation of the
                                                     (momentum) flux),
                         u  = u_D      on domain boundary.

The r.h.s. term g corresponds to gravity, the term f is surface tension
force (here f = kappa · n where kappa is the mean curvature (1/R)). The
Dirichlet data is chosen according to match a manufactured solution
introduced in [1] which allows us to measure errors after the
computation of a discrete solution. The coefficients mu are domain-wise
constants which are different in the two sub-domains.

Discretisation:
---------------
* Finite element space: As in cutfem.py but for every velocity component
  and the pressure.

* Variational formulation: We use a Nitsche formulation which involves
  averages of the fluxes and jumps of the solution across the interface
  [2]. For the average we use the Hansbo-choice [3] where the average
  is adjusted to the local cut configuration in order to ensure stability
  of the resulting viscosity formulation. To ensure inf-sup for the
  velocity-pressure space we add a ghost penalty stabilization on the
  pressure space, cf. [2]. 

* Surface tension: In this example we prescribe the surface tension
  analytically, i.e. we circumvent approximating the mean curvature from
  the level set function. 

Implementational aspects:
-------------------------
* Geometry approximation: As in nxfem_higher_order.py

* Ghost penalty stabilization: The edge-based stabilizations require
  different couplings (compared to standard discretisations). To add
  these to the sparsity pattern we have to add the "dgjumps" flags which
  prepares the sparse matrix for the corresponding needed couplings.

* Linear systems: A (sparse) direct solver is applied to solve the
  arising linear systems.

Literature:
-----------
[1] M. Kirchhart, S. Groß, A. Reusken, Analysis of an XFEM
    discretization for Stokes interface problems, SIAM J. Sci. Comput.,
    38(2), A1019–A1043, 2016.
[2] P. Lederer, C.-M. Pfeiler, C. Wintersteiger, C. Lehrenfeld, Higher
    order unfitted fem for Stokes interface problems. Proc. Appl. Math.
    Mech., 16: 7-10, 2016.
[3] A. Hansbo, P. Hansbo, An unfitted finite element method, based on
    Nitsche's method, for elliptic interface problems, Comp. Meth. Appl.
    Mech. Eng., 191(47):5537-5552, 2002

"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.lsetcurv import *


# -------------------------------- PARAMETERS ---------------------------------
# ----------------------------------- MAIN ------------------------------------

### generate the background mesh
square = SplineGeometry()
square.AddRectangle((-1, -1), (1, 1), bcs=[1, 2, 3, 4])
mesh = Mesh (square.GenerateMesh(maxh=4, quad_dominated=False))
for i in range(4):
    mesh.Refine()
d = mesh.dim


### define problem parameters and data
mu = [1, 10]

R = 2.0 / 3.0

aneg = 1.0/mu[0]
apos = 1.0/mu[1] + (1.0/mu[0] - 1.0/mu[1]) * exp(x*x + y*y - R*R)


## construct the exact level set function, radius and radius squared:
rsqr = x**2 + y**2
r = sqrt(rsqr)

levelset = r - R


## construct the exact solution:
gammaf = 0.5 # surface tension = pressure jump
vel_exact = [a * CoefficientFunction((-y, x)) * exp(-1*rsqr) for a in [aneg, apos]]
pres_exact = [x**3, x**3 - gammaf]

# some helper expressions to compute the r.h.s. of the Stokes system for a given solution:
def Coef_Grad(v):
    return CoefficientFunction(tuple([v[i].Diff(w) for w in [x, y] for i in [0, 1]]), dims=(d, d)).trans
def Coef_Eps(v):
    return 0.5*(Coef_Grad(v) + Coef_Grad(v).trans)
def Coef_Div(v):
    return CoefficientFunction(tuple([sum([v[i, j].Diff(w) for j, w in [(0, x), (1, y)]]) for i in range(d)]))


## manufactured right-hand side:
coef_g = [Coef_Div(-2 * mu[i] * Coef_Eps(vel_exact[i]) + pres_exact[i] * Id(d)) for i in [0, 1]]


## discretization parameters
order = 3 # =k in P(k)/P(k-1)-Taylor-Hood pair

# stabilization parameters for ghost-penalty
gamma_stab_vel = 0.05  # if set to zero: no GP-stabilization for velocity
gamma_stab_pres = 0.05

# stabilization parameter for Nitsche
lambda_nitsche  = 0.5 * (mu[0] + mu[1]) * 20 * order * order


### discretization

## (higher order) level set machinery
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=10.5, discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lsetp1 = lsetmeshadap.lset_p1
ci = CutInfo(mesh, lsetp1)

# finite element space
Vhbase = VectorH1(mesh, order=order, dirichlet=[1, 2, 3, 4])
Qhbase = H1(mesh, order=order-1)

Vhneg = Compress(Vhbase, GetDofsOfElements(Vhbase, ci.GetElementsOfType(HASNEG)))
Vhpos = Compress(Vhbase, GetDofsOfElements(Vhbase, ci.GetElementsOfType(HASPOS)))
Qhneg = Compress(Qhbase, GetDofsOfElements(Qhbase, ci.GetElementsOfType(HASNEG)))
Qhpos = Compress(Qhbase, GetDofsOfElements(Qhbase, ci.GetElementsOfType(HASPOS)))

WhG = FESpace([Vhneg*Vhpos,Qhneg*Qhpos,NumberSpace(mesh)], dgjumps=True)

# discrete functions
gfup = GridFunction(WhG)
gfu, gfp, gfn = gfup.components
gfu_neg, gfu_pos = gfu.components
gfp_neg, gfp_pos = gfp.components

# describe the integration regions
lset_neg = {"levelset": lsetp1, "domain_type": NEG}
lset_pos = {"levelset": lsetp1, "domain_type": POS}
lset_doms = [lset_neg, lset_pos]
lset_if  = {"levelset": lsetp1, "domain_type": IF }

# element, facet and dof marking w.r.t. boundary approximation with lsetp1:
ba_facets = [GetFacetsWithNeighborTypes(mesh, a=ci.GetElementsOfType(HASNEG), b=ci.GetElementsOfType(IF)),
             GetFacetsWithNeighborTypes(mesh, a=ci.GetElementsOfType(HASPOS), b=ci.GetElementsOfType(IF))]


## construct (Bi-)linear forms
h = specialcf.mesh_size 
n_lset = 1.0/Norm(grad(lsetp1)) * grad(lsetp1)
kappa = [CutRatioGF(ci),1.0-CutRatioGF(ci)]
           
a = BilinearForm(WhG, symmetric=False)
f = LinearForm(WhG)
            
u, p, n = WhG.TrialFunction()
v, q, m = WhG.TestFunction()

# define some helper expressions:
def eps(u):
    return 0.5 * (Grad(u) + Grad(u).trans)
def sigma(i, u, p):
    return - 2 * mu[i] * eps(u[i]) + p[i] * Id(mesh.dim)
def average_flux(u, p):
    return sum([kappa[i] * sigma(i, u, p) * n_lset for i in range(2)])
def average_inv(u):
    return sum([kappa[1-i] * u[i] for i in range(2)])
def jump(u):
    return u[0] - u[1]

# Stokes variational formulation:
for i in [0, 1]: # loop over domains
    a += SymbolicBFI(lset_doms[i], form = 2 * mu[i] * InnerProduct(eps(u[i]), eps(v[i])))
    a += SymbolicBFI(lset_doms[i], form = - div(u[i]) * q[i] - div(v[i]) * p[i])
    f += SymbolicLFI(lset_doms[i], form = coef_g[i] * v[i])
# constraining the pressure integral
a += SymbolicBFI(lset_neg, form = n * q[0] + m * p[0])
# Nitsche parts:
a += SymbolicBFI(lset_if,  form = average_flux(u,p) * jump(v))
a += SymbolicBFI(lset_if,  form = average_flux(v,p) * jump(u))
a += SymbolicBFI(lset_if,  form = lambda_nitsche/h * jump(u) * jump(v))
# surface tension:
f += SymbolicLFI(lset_if,  form = -gammaf * average_inv(v)*n_lset)

# ghost-penalty terms:
for i in [0, 1]:
    if gamma_stab_vel > 0:
        a += SymbolicFacetPatchBFI(form = gamma_stab_vel/h**2 * (u[i] - u[i].Other()) * (v[i] - v[i].Other()),
                                   skeleton=False, definedonelements=ba_facets[i])
    if gamma_stab_pres > 0:
        a += SymbolicFacetPatchBFI(form = -gamma_stab_pres * (p[i] - p[i].Other()) * (q[i] - q[i].Other()),
                                   skeleton=False, definedonelements=ba_facets[i])


## apply mesh adaptation for higher order geometry approximation
mesh.SetDeformation(deformation)


### assemble and solve the resulting linear system
with TaskManager():
    a.Assemble()
    f.Assemble()

    # set the in-homogeneous Dirichlet condition
    gfu.components[1].Set(vel_exact[1])

    # solve linear system
    f.vec.data -= a.mat * gfup.vec
    
    gfup.vec.data += a.mat.Inverse(WhG.FreeDofs()) * f.vec
              
# measure the error
vl2error = sqrt(sum([Integrate(lset_doms[i], Norm(gfu.components[i] - vel_exact[i])**2, mesh=mesh) for i in [0, 1]]))
vh1error = sqrt(sum([Integrate(lset_doms[i], Norm(Grad(gfu.components[i]) - Coef_Grad(vel_exact[i]))**2, mesh=mesh) for i in [0, 1]]))
pl2error = sqrt(sum([Integrate(lset_doms[i], Norm(gfp.components[i] - pres_exact[i])**2, mesh=mesh) for i in [0, 1]]))

print("L2 Error of velocity: {:10.8e}".format(vl2error))
print("H1 Error of velocity: {:10.8e}".format(vh1error))
print("L2 Error of pressure: {:10.8e}".format(pl2error))

# unset mesh adaptation
mesh.UnsetDeformation()

# visualization of the solution
Draw(deformation, mesh, "deformation")
Draw(levelset, mesh, "levelset")
Draw(lsetp1, mesh, "lsetp1")
Draw(IfPos(lsetp1, gfu.components[1], gfu.components[0]), mesh, "vel")
Draw(IfPos(lsetp1, gfp.components[1], gfp.components[0]), mesh, "pressure")
