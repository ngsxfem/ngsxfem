"""
In this example we solve a 2-dimensional scalar-valued Laplace-Beltrami
problem on a surface (unit circle) embedded in a square domain. As to
discretisation method we use a levelset-based implicit geometry
description, and a trace finite element method with a normal diffusion
for consistent stabilization. To compute the error we construct a
manufactured solution of the PDE to study optimal convergence for
different polynomial degrees.

Domain:
-------

The background domain Omega is [-1.5,1.5]^2 in R^2 while the surface
Gamma is closed, fixed in space and described implicitly via the zero
level of a level set function (unit circle). For the discretisation the
level set function is approximated first-order by a piecewise linear
interpolation, while the mesh deformation is applied to the geometry for
high-order accuracy.

PDE problem:
------------
Domain equation:
                 u - div_Gamma(grad_Gamma(u)) = f    on Gamma

where div_Gamma(grad_Gamma(( )) is the Laplace-Beltrami operator defined
on a closed surface Gamma which is embedded in a domain Omega in R^d.
The r.h.s. term f is chosen f = 2*(x + y) according to a manufactured
solution u = x + y which allows us to measure errors after the
computation of a discrete solution.

Discretisation:
---------------
* Finite element space: We consider a given, fixed mesh on Omega and a
  standard finite element space Vh on Omega. For a start we use the
  space of piecewise linear functions for Vh. To obtain a finite element
  space for the PDE problem on Gamma we use the restriction of Vh on
  Gamma, the trace finite element Vh^Gamma = tr_Gamma Vh.

* Variational formulation: We use normal diffusion method for consistent
  stabilization in terms of the bilinear form of grad(u)*n and grad(v)*n,
  in order to ensure stability of the resulting formulation.

* Linear systems: A (sparse) direct solver is applied to solve the
  arising linear systems.

Extensions:
-----------
* Instead of using the normal diffusion method, one could use ghost
  penalty method (or some others) for consistent stabilization.
* We choose the order of mesh deformation the same as the finite element
  polynomial degree, which gives the optimal error convergence. You can
  also set the order different from the polynomial degree, in order to
  test other possibilities.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from netgen.csg import *
from ngsolve import *
from ngsolve.internal import *
from xfem import *
from xfem.lsetcurv import *

import matplotlib.pyplot as plt


# -------------------------------- PARAMETERS ---------------------------------
# Domain corners
ll, ur = (-1.5, -1.5), (1.5, 1.5)
# mesh size
maxh = 0.5
# Finite element order
order = 2

# times of mesh refinement for convergence study
n_refs = 5


# ----------------------------------- MAIN ------------------------------------
# a list to store errors
errors = []

# We generate the background mesh of the domain and use a simplicial
# triangulation to obtain a mesh with quadrilaterals use
# 'quad_dominated=True'
square = SplineGeometry()
square.AddRectangle(ll, ur, bc=1)
mesh = Mesh(square.GenerateMesh(maxh=maxh, quad=False))
Draw(mesh)

# Generate a level set function of unit circle
phi = sqrt(x**2 + y**2) - 1.0
Draw(phi, mesh, "levelset")

# Iteration with mesh refinement (bisection) for convergence study
for i in range(n_refs):
    if i > 0:
        mesh.Refine()

    Draw(phi, mesh, "levelset")

    # level set function: 1st-order interpolation and
    #                    'order'th-order deformation
    # threshold: the limit of deformation
    lsetad = LevelSetMeshAdaptation(mesh, order=order, threshold=1000)
    deform = lsetad.CalcDeformation(phi)
    lsetip = lsetad.lset_p1  # InterpolateToP1(lset,lsetip)
    subdiv = 0
    Draw(lsetip, mesh, "lsetapprox")

    # Set mesh deformation to reach geometrically high order accuracy
    mesh.SetDeformation(deform)

    # Declare the integration domains
    lsetif = {"levelset": lsetip, "domain_type": IF, "subdivlvl": subdiv}

    # Declare a trace finite element space with polynomial degree p
    trVh = H1(mesh, order=order, dirichlet=[])
    # Declare the symbolic trial and test functions
    u = trVh.TrialFunction()
    v = trVh.TestFunction()

    # Collect the information about the cut elements
    cut = CutInfo(mesh, lsetip)
    elem = cut.GetElementsOfType(IF)
    cutdof = trVh.FreeDofs() & GetDofsOfElements(trVh, elem)

    # Declare a grid function to store the solution
    gf = GridFunction(trVh)

    # Calculate normal vector n and mesh size h
    n = 1.0 / sqrt(InnerProduct(grad(lsetip), grad(lsetip))) * grad(lsetip)
    h = specialcf.mesh_size
    gamma = 1.0 / h

    # Define the tangential projection P
    def P(u):
        return u - (u * n) * n

    # Assemble the bilinear form A(u,v)
    a = BilinearForm(trVh, symmetric=True, check_unused=False)
    a += SymbolicBFI(levelset_domain=lsetif,
                     form=u * v + P(grad(u)) * P(grad(v)),
                     definedonelements=elem)
    a += SymbolicBFI(form=gamma * (grad(u) * n) * (grad(v) * n),
                     definedonelements=elem)  # normal diffusion
    a.Assemble()

    # Assemble the linear form f(v)
    f = LinearForm(trVh)
    f += SymbolicLFI(levelset_domain=lsetif, form=2 * (x + y) * v,
                     definedonelements=elem)
    f.Assemble()

    # Solve the linear system Au = f
    gf.vec[:] = 0.0
    gf.vec.data = a.mat.Inverse(cutdof) * f.vec
    num = gf

    # Visualization settings
    visoptions.mminval = 0
    visoptions.mmaxval = 0
    visoptions.deformation = 0
    visoptions.autoscale = 1

    # Compute the error between the numerical and the exact solutions
    exa = CoefficientFunction((x + y))
    error = sqrt(Integrate(lsetif, (num - exa) * (num - exa), mesh))
    print("L2 Error: {:.5e}".format(error))

    RefineAtLevelSet(lsetad.lset_p1)
    errors.append((i, error))

    Draw(exa, mesh, "exact")
    Draw(num, mesh, "numeric")

    mesh.UnsetDeformation()
    input("<press enter to continue>")

# plot errors to show convergence
ref = [[errors[0][1] * 0.5**(j * L) for L in range(len(errors))]
       for j in range(n_refs)]

plt.yscale('log')
plt.xlabel("lvl")
plt.ylabel("L2 error")
lvl, err = zip(*errors)
plt.plot(lvl, err, "-*")

for j in range(n_refs):
    plt.plot(lvl, ref[j])
plt.title("order = " + str(order))
plt.ion()
plt.show()
plt.savefig("order = " + str(order) + ".png")
input("<press enter to quit>")
