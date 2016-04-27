# integration on lset domains
from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

from netgen.csg import CSGeometry, OrthoBrick, Pnt
from netgen.meshing import MeshingParameters


# cheese sphere geometry:
# 'Dziuk, Elliott, Finite element methods for surface PDEs, Acta Numerica, 2013', pp. 373-374:
cube = CSGeometry()
cube.Add (OrthoBrick(Pnt(-2.5,-2.5,-2.5), Pnt(2.5,2.5,2.5)))
mesh = Mesh(cube.GenerateMesh (maxh=0.8))
levelset = (sqrt((x*x-1)*(x*x-1)+(y*y-1)*(y*y-1)+(z*z-1)*(z*z-1)+(x*x+y*y-4)*(x*x+y*y-4)+(x*x+z*z-4)*(x*x+z*z-4)+(y*y+z*z-4)*(y*y+z*z-4))-4).Compile()

# dziuk elliott geometry:
# 'Dziuk, Elliott, Finite element methods for surface PDEs, Acta Numerica, 2013', pp. 318-319:
# cube = CSGeometry()
# cube.Add (OrthoBrick(Pnt(-2.5,-1.5,-1.5), Pnt(2.5,1.5,1.5)))
# mesh = Mesh(cube.GenerateMesh (maxh=0.8))
# levelset = (sqrt(0.25*x*x+y*y+4.0*z*z/((1+0.5*sin(pi*x))*(1+0.5*sin(pi*x))))-1.0).Compile()

order = 2
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=10.2, discontinuous_qn=True)

Draw(levelset,mesh,"levelset")
Draw(lsetmeshadap.deform,mesh,"deformation")
Draw(lsetmeshadap.lset_p1,mesh,"levelset(P1)")

distances = []
for reflevel in range(3):

    if(reflevel > 0):
        mesh.Refine()

    # Applying the mesh deformation
    deformation = lsetmeshadap.CalcDeformation(levelset)

    mesh.SetDeformation(deformation)
    distances.append(lsetmeshadap.CalcMaxDistance(levelset));
    mesh.UnsetDeformation()

    # refine cut elements:
    RefineAtLevelSet(gf=lsetmeshadap.lset_p1)

    eoc = [ log(distances[i-1]/distances[i])/log(2) for i in range(1,len(distances))]
    print("distances = {}".format(distances))
    print("eoc = {}".format(eoc))
    # input("press ENTER to refine")
