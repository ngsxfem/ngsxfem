from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube

def test_cutinfo():
    mesh = Mesh("pytests/mesh2D.vol.gz")

# mesh = Mesh("mesh2D.vol.gz")
# lset = GridFunction(H1(mesh,order=1))
# lset.Set(sqrt(x*x+y*y)-0.8)
# ci = CutInfo(mesh)
# ci.Update(lset)
# ba_vol = [ci.GetElementsOfType(dt,VOL) for dt in [NEG,POS,IF]]
# ba_bnd = [ci.GetElementsOfType(dt,BND) for dt in [NEG,POS,IF]]
# for dt in [NEG,POS,IF]:
#     print(ci.GetElementsOfType(dt,VOL))
# for dt in [NEG,POS,IF]:
#     print(ci.GetElementsOfType(dt,BND))
# print(ci.GetCutRatios(VOL))
# print(ci.GetCutRatios(BND))

# Draw(lset,mesh,"lset")
