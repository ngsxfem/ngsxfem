from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.la import *

from libngsxfem_py.xfem import *
from deform import *

from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

geom = SplineGeometry("d01_testgeom.in2d")
mp = MeshingParameters (maxh=0.1)
mesh = Mesh(geom.GenerateMesh (mp))

# lset_stats = StatisticContainer()
# deform_stats = StatisticContainer()
order = 4

R = 2.0/3.0
levelset = VariableCF("sqrt(x*x+y*y)-(0.5+0.1*sin(8*atan2(x,y)))")

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2)

error_tab = []
error_tab_order_k = [[] for k in range(order+2)]

adaptive = True
moreadaptive = False

nrefs = 35
for cnt in range(nrefs+1):
    if cnt > 0 :
        maxdist = lsetmeshadap.CalcMaxDistance(levelset);
        error_tab.append ( maxdist )
        if (cnt==1):
            c = [200 * maxdist / pow(2,k) for k in range(order+2)]
        for k in range(order+2):
            error_tab_order_k[k].append ( c[k] / pow(pow(2.0,k),cnt) )
        print(maxdist)
        if (adaptive):
            if moreadaptive:
                lsetmeshadap.MarkForRefinement(levelset,0.25*maxdist,absolute=True)
            else:
                lsetmeshadap.MarkForRefinement(levelset,1e-99,absolute=False)
    if cnt > 0 and maxdist < 1e-12:
        break
    if cnt != nrefs - 1:
        mesh.Refine()
        deformation = lsetmeshadap.CalcDeformation(levelset)
Draw(deformation)
Draw(lsetmeshadap.lset_p1)

import matplotlib.pyplot as plt
plt.yscale('log')
plt.plot(error_tab, "-*", label="max dist (k = " + str(order) + ")")
for k in range(max(order-1,1),order+2):
    plt.plot(error_tab_order_k[k], "-+", label="order "+str(k))
plt.legend()
plt.ion()
plt.show()
 

input("finished")  
