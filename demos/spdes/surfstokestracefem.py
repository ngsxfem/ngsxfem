"""
In this example we solve a scalar Laplace-Beltrami problem with a
similar discretisation method to the one used in tracefem.py. However,
we use a 3D (background mesh dimension) problem and higher order method
this time. To be robust w.r.t. the interface position also in the
condition number we use the normal diffusion stabilization, cf.[1,2].

Used features:
--------------
* Higher order geometry approximation, cf. jupyter tutorial `lsetint`

* Restricted finite element space to condense the system to active dofs,
  cf. jupyter tutorial `basics`

* Visualization: The visualization of the solution is most convenient
  with paraview and the generated vtk file.

Literature:
-----------
[1] J. Grande, C. Lehrenfeld, A. Reusken, Analysis of a high-order trace
    finite element method for PDEs on level set surfaces.
    SIAM Journal on Numerical Analysis, 56(1):228-255, 2018.

"""

# ------------------------------ LOAD LIBRARIES -------------------------------
import sys

from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve import *
from ngsolve.internal import *
from xfem import *
from xfem.lsetcurv import *
from math import pi

# -------------------------------- PARAMETERS ---------------------------------
# Mesh diameter
maxh = 0.6
# Geometry
geom = "circle"
# Polynomial order of FE space
order = 2

# Problem parameters
reac_cf = 1
diff_cf = 1
# Ellipsoid parameter
c = 1
# Geometry
cube = CSGeometry()
if geom == "circle":
    phi = Norm(CoefficientFunction((x / c, y, z))) - 1

    cube.Add(OrthoBrick(Pnt(-2, -2, -2), Pnt(2, 2, 2)))
    mesh = Mesh(cube.GenerateMesh(maxh=maxh, quad_dominated=False))
elif geom == "decocube":
    phi = (x ** 2 + y ** 2 - 4) ** 2 + (y ** 2 - 1) ** 2 + (y ** 2 + z ** 2 - 4) ** 2 + (x ** 2 - 1) ** 2 + (
                x ** 2 + z ** 2 - 4) ** 2 + (z ** 2 - 1) ** 2 - 13
    cube.Add(OrthoBrick(Pnt(-3, -3, -3), Pnt(3, 3, 3)))
    mesh = Mesh(cube.GenerateMesh(maxh=maxh, quad_dominated=False))

# Preliminary refinements. Keep it small as assembling the linear form is costly if the exact solution is not tangential
n_cut_ref = 0
for i in range(n_cut_ref):
    lsetp1 = GridFunction(H1(mesh, order=2))
    InterpolateToP1(phi, lsetp1)
    RefineAtLevelSet(lsetp1)
    mesh.Refine()
# Class to compute the mesh transformation needed for higher order accuracy
#  * order: order of the mesh deformation function
#  * threshold: barrier for maximum deformation (to ensure shape regularity)

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(phi)
lset_approx = lsetmeshadap.lset_p1

mesh.SetDeformation(deformation)
# Class to compute the mesh transformation needed for higher order accuracy
#  * order: order of the mesh deformation function
#  * threshold: barrier for maximum deformation (to ensure shape regularity)

# Background FESpaces
Vh = VectorH1(mesh, order=order, dirichlet=[])
Qh = H1(mesh, order=order - 1, dirichlet=[])
ci = CutInfo(mesh, lset_approx)
ba_IF = ci.GetElementsOfType(IF)
VhG = Restrict(Vh, ba_IF)
L2G = Restrict(Qh, ba_IF)

X = VhG * L2G
gfu = GridFunction(X)


# Helper Functions for tangential projection and gradients of the exact solution
def coef_grad(u):
    dirs = {1: [x], 2: [x, y], 3: [x, y, z]}
    if u.dim == 1:
        return CF(tuple([u.Diff(r) for r in dirs[mesh.dim]]))
    else:
        return CF(tuple([u[i].Diff(r) for i in range(u.dim) for r in dirs[mesh.dim]]), dims=(u.dim, mesh.dim))


def grad(u):
    if type(u) in [ProxyFunction, GridFunction]:
        return ngsolve.grad(u)
    else:
        return coef_grad(u)


def Pmats(n):
    return Id(3) - OuterProduct(n, n)


# Coefficients / parameters:

n = Normalize(grad(lset_approx))
ntilde = Interpolate(grad(phi), VhG)

n.Compile()
ntilde.Compile()
h = specialcf.mesh_size

eta = 100.0 / (h * h)
rhou = 1.0 / h
rhop = h
# Exact tangential projection
Ps = Pmats(Normalize(grad(phi)))
# define solution and right-hand side


if geom == "circle":
    functions = {
        "extvsol1": ((-y - z) * x + y * y + z * z),
        "extvsol2": ((-x - z) * y + x * x + z * z),
        "extvsol3": ((-x - y) * z + x * x + y * y),
        "rhs1": -((y + z) * x - y * y - z * z) * (x * x + y * y + z * z + 1) / (x * x + y * y + z * z),
        "rhs2": ((-x - z) * y + x * x + z * z) * (x * x + y * y + z * z + 1) / (x * x + y * y + z * z),
        "rhs3": ((-x - y) * z + x * x + y * y) * (x * x + y * y + z * z + 1) / (x * x + y * y + z * z),
    }
    extpsol = (x * y ** 3 + z * (x ** 2 + y ** 2 + z ** 2) ** (3 / 2)) / ((x ** 2 + y ** 2 + z ** 2) ** 2)
    uSol = Ps * CoefficientFunction((functions["extvsol1"], functions["extvsol2"], functions["extvsol3"]))
    pSol = CoefficientFunction((functions["extvsol1"], functions["extvsol2"], functions["extvsol3"]))

elif geom == "decocube":
    uSol = Ps * CoefficientFunction((-z ** 2, y, x))
    pSol = CoefficientFunction(
        x * y ** 3 + z - 1 / Integrate({"levelset": lset_approx, "domain_type": IF}, cf=CoefficientFunction(1),
                                       mesh=mesh) * Integrate({"levelset": lset_approx, "domain_type": IF},
                                                              cf=x * y ** 3 + z, mesh=mesh))
    extpsol = pSol
    extpsol.Compile()

u, p = X.TrialFunction()
v, q = X.TestFunction()
# Measure on surface
ds = dCut(lset_approx, IF, definedonelements=ba_IF, deformation=deformation)
# Measure on the bulk around the surface
dx = dx(definedonelements=ba_IF, deformation=deformation)


# Helper Functions for rate-of-strain-tensor
def eps(u):
    return Ps * Sym(grad(u)) * Ps


def divG(u):
    if u.dim == 3:
        return Trace(grad(u) * Ps)
    if u.dim == 9:
        N = 3
        divGi = [divG(u[i, :]) for i in range(N)]
        return CF(tuple([divGi[i] for i in range(N)]))


# Weingarten mappings
weing = grad(n)
weingex = grad(Normalize(grad(phi)))

with TaskManager():
    rhs1 = Ps * -divG(eps(uSol) - (uSol * Normalize(grad(phi))) * weingex) + uSol + grad(extpsol)

    rhs2 = Trace(Ps * grad(uSol))

    # bilinear forms:
    Pmat = Pmats(n)


    def E(u):
        return Pmat * Sym(grad(u)) * Pmat - (u * n) * weing


    a = BilinearForm(X, symmetric=True)
    a += InnerProduct(E(u), E(v)) * ds
    a += (Pmat * u) * (Pmat * v) * ds
    a += (eta * ((u * ntilde) * (v * ntilde))) * ds
    a += rhou * ((grad(u) * n) * (grad(v) * n)) * dx
    a += InnerProduct(u, Pmat * grad(q)) * ds
    a += InnerProduct(v, Pmat * grad(p)) * ds
    a += -rhop * InnerProduct(n * grad(p), n * grad(q)) * dx

    a.Assemble()

    # R.h.s. linear form
    f = LinearForm(X)
    f += rhs1 * v * ds
    f += -rhs2 * q * ds
    f.Assemble()

    gfu.vec.data = a.mat.Inverse(freedofs=X.FreeDofs()) * f.vec

    print("l2error : ", sqrt(Integrate((gfu.components[0] - uSol) ** 2 * ds, mesh=mesh)))
uerr = (gfu.components[0] - uSol)

Draw(deformation, mesh, "deformation")
Draw(gfu.components[0], mesh, "u")

if len(sys.argv) < 2 or not (sys.argv[1] == "testmode"):
    input("Continue (press enter) to create a VTK-Output to stokestracefem3d.vtk")

    with TaskManager():
        vtk = VTKOutput(ma=mesh,
                        coefs=[lset_approx, gfu.components[0], uSol, uerr],
                        names=["P1-levelset", "uh", "uex", "uerr"],
                        filename="stokestracefem3d", subdivision=2)

        vtk.Do()
