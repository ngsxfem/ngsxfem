#from ngsolve.fem import *
#from ngsolve.comp import *
#from ngsolve.solve import *
#from ngsolve.la import *
#from ngsolve.utils import *
#
#from xfem import *

from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *



from ctypes import CDLL
# on Windows replace '.so' with '.dll'
mylngs = CDLL("libngsxfem_spacetime.so")


square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
ngmesh = square.GenerateMesh(maxh=0.25, quad_dominated=False)
mesh = Mesh (ngmesh)

fes1 = V=H1(mesh, order=1, dirichlet=[1,2,3,4])
tfe = ScalarTimeFE(1);

st_fe = SpaceTimeFESpace(fes1,tfe)
st_fe.SetTime(0.5)
#fes = FESpace("myfespace", mesh, order = 1, dirichlet=[1,2,3,4], flags = {"order_time" : True } )

#help(fes1)
#input("")
#print ("freedofs: ", fes.FreeDofs())

visoptions.autoscale = False
visoptions.mminval=0.0
visoptions.mmaxval=1.0
visoptions.deformation = 1

w = GridFunction(st_fe)

#L = [el for el in fes.Elements()]
#print("Len(L) = {0}".format(len(L)))
#print("Dof-Nrs. = {0}".format(L[0].dofs))

Draw(w)

for i in range(st_fe.ndof):
   # print("i = {0}".format(i))
    if i > 0:
        w.vec[i-1] = 0
    w.vec[i] = 1
    Redraw()
    #input("")
    sleep(0.25)



