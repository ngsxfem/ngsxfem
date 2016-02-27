# solve the Laplace Beltrami equation on a sphere
# with Dirichlet boundary condition u = 0

# ngsolve stuff
from ngsolve import *

# xfem and trace fem stuff
import libngsxfem_py.xfem as xfem                                 

# geometry stuff
# unit square [-2,2]^2 (boundary index 1 everywhere)
def MakeSquare():
    from netgen.geom2d import SplineGeometry
    square = SplineGeometry()
    square.AddRectangle([-2,-2],[2,2],bc=1)
    return square;

# unit cube [-2,2]^3 (boundary index 1 everywhere)
def MakeCube():
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    cube = CSGeometry()
    cube.Add (OrthoBrick(Pnt(-2,-2,-2), Pnt(2,2,2)))
    return cube;

# meshing
from netgen.meshing import MeshingParameters
ngsglobals.msg_level = 1
mesh = Mesh (MakeCube().GenerateMesh(maxh=0.5, quad_dominated=False))

# The mesh deformation calculator
from deform import *
order = 4
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)

# The level set
levelset = sqrt(x*x+y*y) - 1

# Calculation of the deformation:
deformation = lsetmeshadap.CalcDeformation(levelset)
print ("Computed mesh transformation.")
lset_ho = lsetmeshadap.lset_ho
# lsetmeshadap.lset_p1
print ("Max. dist of discr. intf. to exact level set fct:", lsetmeshadap.CalcMaxDistance(levelset))
print ("Max. dist of discr. intf. to discr. level set f.:", lsetmeshadap.CalcMaxDistance(lset_ho))

# Applying the mesh deformation
mesh.SetDeformation(deformation)

Vh = H1(mesh, order=3, dirichlet=[1])
Vh_tr = CastToXFESpace( FESpace ("xfespace", mesh=mesh, order=3, 
                                       flags = {"trace" : True, "dgjumps" : False, "ref_space" : 0}))
Vh_tr.SetBaseFESpace(Vh)
Vh_tr.SetLevelSet(levelset)
Vh_tr.Update()

print ("Trace FE Space created.")
#unset the mesh deformation
mesh.UnsetDeformation() 
mesh.Refine()


print ("ran through.")
