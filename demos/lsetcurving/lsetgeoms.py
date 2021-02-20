"""
This examples shows a number of pre-implemented 3d level set geometries,
computes the distance between the smooth level set and the P1
approximation as well as the order of convergence thereof.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.csg import CSGeometry
from ngsolve import *
from ngsolve.internal import *

from xfem import *
from xfem.lsetcurv import *
from xfem.utils import LevelsetExamples, BoundingBoxes

import sys


# -------------------------------- PARAMETERS ---------------------------------
maxh = 0.5
order = 2
maxreflvl = 4


# ----------------------------------- MAIN ------------------------------------
# Visual options for the Netgen GUI
# viewoptions.clipping.enable = 1
# viewoptions.clipping.nx = 0.0
# viewoptions.clipping.ny = 1
# viewoptions.clipping.nz = 0.0
# visoptions.mminval = 0.0
# visoptions.mmaxval = 0.0
# visoptions.autoscale = False
# visoptions.isosurf = 1
# visoptions.numiso = 1
# visoptions.subdivisions = 1

# Main loop
for lsetgeom in ["cheese", "torus", "dziukelliott", "dziuk88", "sphere"]:
    print('Geometry: ', lsetgeom)

    geom = CSGeometry()
    geom.Add(BoundingBoxes[lsetgeom])
    mesh = Mesh(geom.GenerateMesh(maxh=maxh))
    levelset = LevelsetExamples[lsetgeom]

    lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=100,
                                          discontinuous_qn=True)

    Draw(levelset, mesh, "levelset")
    Draw(lsetmeshadap.deform, mesh, "deformation")
    Draw(lsetmeshadap.lset_p1, mesh, "levelset(P1)")

    distances = []
    for reflevel in range(maxreflvl):

        if(reflevel > 0):
            mesh.Refine()

        deformation = lsetmeshadap.CalcDeformation(levelset)
        distances.append(lsetmeshadap.CalcMaxDistance(levelset, deform=True))
        
        # Refine cut elements:
        RefineAtLevelSet(gf=lsetmeshadap.lset_p1)

        # Post-processing
        Redraw(blocking=True)
        eoc = [log(distances[i - 1] / distances[i]) / log(2)
               for i in range(1, len(distances))]
        print("distances = {}".format(distances))
        print("eoc = {}".format(eoc))

    if (hasattr(sys, 'ps1') or sys.flags.interactive
            or 'running_in_netgen' in dir()):
        input("press ENTER to go to next geometry (or finish)")
