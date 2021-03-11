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
[2] E. Burman, P. Hansbo, M. G. Larson, A. Massing, Cut  finite  element
    methods  for  partial  differential  equations  on  embedded  manifolds
    of  arbitrary  codimensions.  ESAIM: M2AN, 52(6):2247-2282, 2018.
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
# Bisect cut elements of (initial) mesh
n_cut_ref = 2
# Polynomial order of FE space
order = 2

# Problem parameters
reac_cf = 1
diff_cf = 1

# Geometry
cube = CSGeometry()
cube.Add(OrthoBrick(Pnt(-1.41, -1.41, -1.41), Pnt(1.41, 1.41, 1.41)))
mesh = Mesh(cube.GenerateMesh(maxh=maxh, quad_dominated=False))

levelset = sqrt(x**2 + y**2 + z**2) - 1
exact = sin(pi * z)
coef_f = (sin(pi * z) * (diff_cf * pi * pi * (1 - z * z) + reac_cf)
          + diff_cf * cos(pi * z) * 2 * pi * z)

# ----------------------------------- MAIN ------------------------------------
# Preliminary refinements
for i in range(n_cut_ref):
    lsetp1 = GridFunction(H1(mesh, order=1))
    InterpolateToP1(levelset, lsetp1)
    RefineAtLevelSet(lsetp1)
    mesh.Refine()

# Class to compute the mesh transformation needed for higher order accuracy
#  * order: order of the mesh deformation function
#  * threshold: barrier for maximum deformation (to ensure shape regularity)
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000,
                                      discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lset_approx = lsetmeshadap.lset_p1

# Background FESpace
Vh = H1(mesh, order=order, dirichlet=[])

ci = CutInfo(mesh, lset_approx)
ba_IF = ci.GetElementsOfType(IF)
VhG = Restrict(Vh, ba_IF)

gfu = GridFunction(VhG)

# Coefficients / parameters:
n = Normalize(grad(lset_approx))
h = specialcf.mesh_size


# Tangential projection
def P(u):
    return u - (u * n) * n


u, v = VhG.TnT()

# Measure on surface
ds = dCut(lset_approx, IF, definedonelements=ba_IF, deformation=deformation)
# Measure on the bulk around the surface
dx = dx(definedonelements=ba_IF, deformation=deformation)

# Bilinear forms:
a = BilinearForm(VhG, symmetric=True)
a += (diff_cf * P(grad(u)) * P(grad(v)) + reac_cf * u * v) * ds
a += ((diff_cf / h + reac_cf * h) * (grad(u) * n) * (grad(v) * n)) * dx

# R.h.s. linear form
f = LinearForm(VhG)
f += coef_f * v * ds

a.Assemble()
f.Assemble()

gfu.vec.data = a.mat.Inverse(VhG.FreeDofs()) * f.vec

print("l2error : ", sqrt(Integrate((gfu - exact)**2 * ds, mesh=mesh)))
Draw(deformation, mesh, "deformation")
Draw(gfu, mesh, "u")

if len(sys.argv) < 2 or not (sys.argv[1] == "testmode"):
    input("Continue (press enter) to create a VTK-Output to tracefem3d.vtk")
    VTKOutput(ma=mesh, coefs=[deformation, lset_approx, gfu],
              names=["deformation", "P1-levelset", "u"],
              filename="tracefem", subdivision=2).Do()
