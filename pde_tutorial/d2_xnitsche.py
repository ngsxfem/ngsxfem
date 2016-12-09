
# interactive modifications to d2_xnitsche.py
print ("hello from d2_xnitsche.py ;-)")

from ngsolve.solve import *
from ngsolve.comp import *
from ngsolve.fem import *
from ngsolve.la import *
from ngsolve.bla import *
import ngsolve.ngstd as ngstd
from ngsolve.solve import Redraw

from xfem import *

from math import sin                                                  
from time import sleep

def PrintDofs(pde,mesh,fes):
    print ("Printing dofs per element:\n\n")
    for i in mesh.Elements():
        print("dofnrs of element", i, ":\n", fes.GetDofNrs(i))

def ShapeTest(pde,u):
    print ("Shape test:\n")
    u[:][:] = 0                                                   
    for i in range(u.size):                                       
        print ("i = ", i ,".")
        u[:][i-1] = 0.0
        u[:][i] = 1.0                                            
        Redraw(blocking=True)                                      
        sleep(1)

def Test(pde):
    PrintDofs(pde,pde.Mesh(),pde.spaces["fescomp"].StdFESpace)
    PrintDofs(pde,pde.Mesh(),pde.spaces["fescomp"].XFESpace)
    ShapeTest(pde,pde.gridfunctions["u"].vec)
