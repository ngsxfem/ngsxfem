"""
In this example we solve a scalar Laplace-Beltrami problem with a
similar discretisation method to the one used in tracefem.py. However,
we use a 3D (background mesh dimension) problem and higher order method
this time.

Further comments:
-------------------------
* Geometry approximation: We use the same approach as described in
  nxfem_higher_order.py

* Linear systems: To be robust w.r.t. the interface position also in the
  condition number we use the normal diffusion stabilization.

* Visualization: The visualization of the solution is most convenient
  with paraview and the generated vtk file.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve import *
from ngsolve.internal import *
from xfem import *
from xfem.lsetcurv import *

from math import pi


# -------------------------------- PARAMETERS ---------------------------------
# Mesh diameter
maxh = 0.6
# Refine cut elements of (initial) mesh
n_cut_ref = 1
# Polynomial order of FE space
order = 3
# Subdivisions of cut elements in construction of quadrature rule
subdivlvl = 0

# Problem parameters
reac_cf = 1
diff_cf = 1


# ----------------------------------- MAIN ------------------------------------
# Geometry
cube = CSGeometry()
cube.Add(OrthoBrick(Pnt(-1.41, -1.41, -1.41), Pnt(1.41, 1.41, 1.41)))
mesh = Mesh(cube.GenerateMesh(maxh=maxh, quad_dominated=False))

levelset = sqrt(x * x + y * y + z * z) - 1

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

# The piecewise linear interpolation used in the level set adaptation class
lset_approx = lsetmeshadap.lset_p1


# Extended FESpace
VhG = H1(mesh, order=order, dirichlet=[])

# Overwrite freedofs of VhG to mark only dofs that are involved in the
# cut problem
ci = CutInfo(mesh, lset_approx)
ba_IF = ci.GetElementsOfType(IF)
cf_IF = BitArrayCF(ba_IF)
freedofs = VhG.FreeDofs()
freedofs &= GetDofsOfElements(VhG, ba_IF)

gfu = GridFunction(VhG)

# coefficients / parameters:
n = 1.0 / Norm(grad(lset_approx)) * grad(lset_approx)
h = specialcf.mesh_size


# Tangential projection
def P(u):
    return u - (u * n) * n


# Expressions of test and trial functions:
u = VhG.TrialFunction()
v = VhG.TestFunction()

# Integration domains (and integration parameter "subdivlvl" and
# "force_intorder")
lset_if = {"levelset": lset_approx, "domain_type": IF, "subdivlvl": subdivlvl}

# Bilinear forms:
a = BilinearForm(VhG, symmetric=True)
a += SymbolicBFI(levelset_domain=lset_if,
                 form=diff_cf * P(grad(u)) * P(grad(v)) + reac_cf * u * v)
a += SymbolicBFI(
    form=(diff_cf / h + reac_cf * h) * (cf_IF * grad(u) * n) * (grad(v) * n))

# R.h.s. linear form
f_cf = (sin(pi * z) * (diff_cf * pi * pi * (1 - z * z) + reac_cf)
        + diff_cf * cos(pi * z) * 2 * pi * z)

f = LinearForm(VhG)
f += SymbolicLFI(levelset_domain=lset_if, form=f_cf * v)

mesh.deformation = deformation
a.Assemble()
f.Assemble()

gfu.vec[:] = 0.0
gfu.vec.data = a.mat.Inverse(freedofs) * f.vec

exact = sin(pi * z)
err_sqr_coefs = (gfu - exact) * (gfu - exact)
l2error = sqrt(Integrate(levelset_domain=lset_if,
                         cf=err_sqr_coefs[0], mesh=mesh, order=2))

mesh.deformation = None

print("l2error : ", l2error)
Draw(deformation, mesh, "deformation")
Draw(gfu, mesh, "u")

visoptions.mminval = -1
visoptions.mmaxval = 1
visoptions.deformation = 1
visoptions.autoscale = 0

input("Continue (press enter) to create a VTK-Output to tracefem3d.vtk")

vtk = VTKOutput(ma=mesh, coefs=[deformation, lset_approx, gfu],
                names=["deformation", "P1-levelset", "u"],
                filename="tracefem3d", subdivision=2)
vtk.Do()
