from time import sleep
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *

from netgen.geom2d import SplineGeometry
square = SplineGeometry()
square.AddRectangle([-1,-1],[1,1],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.6, quad_dominated=False))

levelset = sqrt(x*x+y*y)-0.7

sleep(1)

# way to create an extended FESpace:
# make the standard space
fes = H1(mesh, order = 1)
# make the extended space
xfes = XFESpace(fes,levelset)
# make a compound from these spaces
xstdfes1 = FESpace([fes,xfes])

u = GridFunction(xstdfes1)

Draw(levelset,mesh,"levelset")
sleep(2)

ci = CutInfo(mesh)
ci.Update(levelset)
Draw(BitArrayCF(ci.GetElementsOfType(IF,VOL)),mesh,"if")
Draw(BitArrayCF(ci.GetElementsOfType(NEG,VOL)),mesh,"neg")
Draw(BitArrayCF(ci.GetElementsOfType(POS,VOL)),mesh,"pos")

Draw(u.components[0]+IfPos(levelset,pos(u.components[1]),neg(u.components[1])),mesh,"u") #sd=6)


for i in range(fes.ndof):
    totali = fes.ndof + i
    u.vec[:][totali-1] = 0.0
    stdi = xfes.BaseDofOfXDof(i)
    u.vec[:][stdi] = 1.0
    print("standard test function {}       ".format(stdi))
    Redraw()
    sleep(2)
    u.vec[:][stdi] = 0.0
    u.vec[:][totali] = 1.0
    print("extended test function {}       ".format(i))
    Redraw()
    sleep(2)
    
