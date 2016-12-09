# integration on lset domains
from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

import sys

# For levelsets
from xfem.utils import LevelsetExamples, BoundingBoxes

from netgen.csg import CSGeometry, OrthoBrick, Pnt
from netgen.meshing import MeshingParameters

# from ngsolve.internal import *
# viewoptions.clipping.enable=1
# viewoptions.clipping.nx = 0.0
# viewoptions.clipping.ny = 1
# viewoptions.clipping.nz = 0.0
# visoptions.mminval=0.0
# visoptions.mmaxval=0.0
# visoptions.autoscale = False
# visoptions.isosurf = 1
# visoptions.numiso = 1
# visoptions.subdivisions = 1

for lsetgeom in ["cheese","torus","dziukelliott","dziuk88","sphere"]:
    geom = CSGeometry()
    geom.Add (BoundingBoxes[lsetgeom])
    mesh = Mesh(geom.GenerateMesh (maxh=1.0))
    levelset = LevelsetExamples[lsetgeom]

    order = 2
    lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=100, discontinuous_qn=True)

    Draw(levelset,mesh,"levelset")
    Draw(lsetmeshadap.deform,mesh,"deformation")
    Draw(lsetmeshadap.lset_p1,mesh,"levelset(P1)")

    distances = []
    for reflevel in range(5):

        if(reflevel > 0):
            mesh.Refine()

        # Applying the mesh deformation
        deformation = lsetmeshadap.CalcDeformation(levelset)
        
        mesh.SetDeformation(deformation)
        distances.append(lsetmeshadap.CalcMaxDistance(levelset));
        mesh.UnsetDeformation()

        # refine cut elements:
        RefineAtLevelSet(gf=lsetmeshadap.lset_p1)

        Redraw(blocking=True)
        eoc = [ log(distances[i-1]/distances[i])/log(2) for i in range(1,len(distances))]
        print("distances = {}".format(distances))
        print("eoc = {}".format(eoc))
    if hasattr(sys,'ps1') or sys.flags.interactive or 'running_in_netgen' in dir():
        input("press ENTER to go to next geometry (or finish)")
