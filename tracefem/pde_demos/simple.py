# solve the Laplace Beltrami equation on a sphere
# with Dirichlet boundary condition u = 0

# ngsolve stuff
from ngsolve import *

# xfem and trace fem stuff
import libngsxfem_py.xfem as xfem                                 

# geometry stuff
# unit square [-1,1]^2 (boundary index 1 everywhere)
# from netgen.geom2d import SplineGeometry
# square = SplineGeometry()
# square.AddRectangle([-1,-1],[1,1],bc=1)

# unit cube [-1,1]^3 (boundary index 1 everywhere)
from netgen.csg import CSGeometry, OrthoBrick, Pnt
cube = CSGeometry()
cube.Add (OrthoBrick(Pnt(-2,-2,-2), Pnt(2,2,2)))

# meshing
from netgen.meshing import MeshingParameters
ngsglobals.msg_level = 1
mesh = Mesh (cube.GenerateMesh(maxh=0.5, quad_dominated=False))

# The mesh deformation calculator
from deform import *
order = 4
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)

# The level set
levelset = sqrt(x*x+y*y) - 1

# Calculation of the deformation:
deformation = lsetmeshadap.CalcDeformation(levelset)
lset_ho = lsetmeshadap.lset_ho
# lsetmeshadap.lset_p1
print ("Max. dist of discr. intf. to exact level set fct:", lsetmeshadap.CalcMaxDistance(levelset))
print ("Max. dist of discr. intf. to discr. level set f.:", lsetmeshadap.CalcMaxDistance(lset_ho))

# Applying the mesh deformation
mesh.SetDeformation(deformation)


#unset the mesh deformation
mesh.UnsetDeformation() 
mesh.Refine()


print ("ran through")
