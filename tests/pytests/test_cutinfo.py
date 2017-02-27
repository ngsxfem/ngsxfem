from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube

def test_cutinfo():
    mesh = Mesh("pytests/mesh2D.vol.gz")

mesh = Mesh("mesh2D.vol.gz")
lset = GridFunction(H1(mesh,order=1))
lset.Set(sqrt(x*x+y*y)-1.0/3.0)
ci = CutInfo(mesh)
ci.Update(lset)
cneg = ci.GetElementsOfType(NEG)
cpos = ci.GetElementsOfType(POS)
cif = ci.GetElementsOfType(IF)
