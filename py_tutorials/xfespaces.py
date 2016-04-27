from time import sleep
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *

from netgen.geom2d import SplineGeometry
square = SplineGeometry()
square.AddRectangle([-1,-1],[1,1],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.6, quad_dominated=False))

levelset = sqrt(x*x+y*y)-0.7

sleep(1)

# ways to create an extended FESpace:
# 1. Way: create standard FESpace and build an XFESpace based on this 

# make the standard space
fes = H1(mesh, order = 1)
# make the extended space
xfes = XFESpace(fes,levelset)
# make a compound from these spaces
xstdfes1 = FESpace([fes,xfes])

# 2. Way: Create an extended standard FESpace

xstdfes2 = XStdFESpace(mesh, levelset, order=1)
u = GridFunction(xstdfes2)

Draw(levelset,mesh,"levelset")
sleep(2)

Draw(u,sd=6)

for i in range(xstdfes2.XFESpace.ndof):
    totali = xstdfes2.StdFESpace.ndof + i
    u.vec[:][totali-1] = 0.0
    stdi = xstdfes2.XFESpace.BaseDofOfXDof(i)
    u.vec[:][stdi] = 1.0
    print("standard test function {}       ".format(stdi))
    Redraw()
    sleep(2)
    u.vec[:][stdi] = 0.0
    u.vec[:][totali] = 1.0
    print("extended test function {}       ".format(i))
    Redraw()
    sleep(2)
    
