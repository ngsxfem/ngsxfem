"""
In this example we solve a scalar Laplace-Beltrami problem with a similar discretization method
to the one used in tracefem.py. However, we use a 3D (background mesh dimension) problem and 
higher order method this time.

    Further comments:
    -------------------------
    Geometry approximation:
    We use the same approach as described in nxfem_higher_order.py

    Linear systems:
    To be robust w.r.t. the interface position also in the condition number we use the
    normal diffusion stabilization.

    Visualization:
    The visualization of the solution is most convenient with paraview and the generated
    vtk file.

"""
from math import pi
# ngsolve stuff
from ngsolve import *
# visualization stuff
from ngsolve.internal import *
# basic xfem functionality
from xfem import *
from xfem.lsetcurv import *

# from netgen.geom2d import SplineGeometry
from netgen.csg import CSGeometry, OrthoBrick, Pnt

# geometry
cube = CSGeometry()
cube.Add (OrthoBrick(Pnt(-1.41,-1.41,-1.41), Pnt(1.41,1.41,1.41)))
mesh = Mesh (cube.GenerateMesh(maxh=0.6, quad_dominated=False))

levelset = sqrt(x*x+y*y+z*z)-1

for i in range(1):
   lsetp1 = GridFunction(H1(mesh,order=1))
   InterpolateToP1(levelset,lsetp1)
   RefineAtLevelSet(lsetp1)
   mesh.Refine()
   
order = 3

# class to compute the mesh transformation needed for higher order accuracy
#  * order: order of the mesh deformation function
#  * threshold: barrier for maximum deformation (to ensure shape regularity)
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)

# the piecewise linear interpolation used in the level set adaptation classe
lset_approx = lsetmeshadap.lset_p1
subdivlvl = 0

# extended FESpace 
VhG = H1(mesh, order=order, dirichlet=[])

# overwrite freedofs of VhG to mark only dofs that are involved in the cut problem
ci = CutInfo(mesh, lset_approx)
ba_IF = ci.GetElementsOfType(IF)
cf_IF = BitArrayCF(ba_IF)
freedofs = VhG.FreeDofs()
freedofs &= GetDofsOfElements(VhG,ba_IF)

gfu = GridFunction(VhG)

# coefficients / parameters: 
n = 1.0/sqrt(InnerProduct(grad(lset_approx),grad(lset_approx))) * grad(lset_approx)
h = specialcf.mesh_size

#tangential projection
def P(u):
   return u - (u*n)*n

# expressions of test and trial functions:
u = VhG.TrialFunction()
v = VhG.TestFunction()

# integration domains (and integration parameter "subdivlvl" and "force_intorder")
lset_if  = { "levelset" : lset_approx, "domain_type" : IF , "subdivlvl" : subdivlvl}

# bilinear forms:
reac_coeff=1
diff_coeff=1
a = BilinearForm(VhG, symmetric = True)
a += SymbolicBFI(levelset_domain = lset_if , form = diff_coeff * P(grad(u)) * P(grad(v)) + reac_coeff * u * v)
a += SymbolicBFI(form = (diff_coeff/h+reac_coeff*h)*(cf_IF * grad(u)*n) * (grad(v)*n))

f_coeff = (sin(pi*z)*(diff_coeff*pi*pi*(1-z*z)+reac_coeff)+diff_coeff*cos(pi*z)*2*pi*z)

f = LinearForm(VhG)
f += SymbolicLFI(levelset_domain = lset_if, form = f_coeff * v)

mesh.SetDeformation(deformation)
a.Assemble()
f.Assemble();

gfu.vec[:] = 0.0
gfu.vec.data = a.mat.Inverse(freedofs) * f.vec

exact = sin(pi*z)
err_sqr_coefs = (gfu-exact)*(gfu-exact)
l2error = sqrt( Integrate( levelset_domain=lset_if, cf=err_sqr_coefs[0], mesh=mesh, order=2) )

mesh.UnsetDeformation()

print ("l2error : ", l2error)
import sys
if not hasattr(sys, 'argv') or len(sys.argv) == 1 or sys.argv[1] != "testmode":
   Draw(deformation,mesh,"deformation")
   Draw(gfu,mesh,"u")

   visoptions.mminval = -1
   visoptions.mmaxval = 1
   visoptions.deformation = 1
   visoptions.autoscale = 0

   input("Continue (press enter) to create a VTK-Output to tracefem3d.vtk")
   
   vtk = VTKOutput(ma=mesh,coefs=[deformation,lset_approx,gfu],names=["deformation","P1-levelset","u"],filename="tracefem3d",subdivision=2)
   vtk.Do()


