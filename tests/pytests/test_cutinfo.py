from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from netgen.geom2d import SplineGeometry

def test_cutinfo():
    mesh = Mesh("pytests/mesh2D.vol.gz")

mesh = Mesh("mesh2D.vol.gz")

# geom = SplineGeometry()
# pnts = [ (0,0), (1,0), (0,1) ]
# pnums = [geom.AppendPoint(*p) for p in pnts]
# lines = [(0,1,1,1,0), (1,2,1,1,0), (2,0,1,1,0)]
# for p1,p2,bc,left,right in lines:
#     geom.Append( ["line", pnums[p1], pnums[p2]], bc=bc, leftdomain=left, rightdomain=right)
# mesh = Mesh (geom.GenerateMesh(maxh=0.05)) # <- will have 4 elements

lset = GridFunction(H1(mesh,order=1))
lset.Set(sqrt(x*x+y*y)-0.8)
ci = CutInfo(mesh)
ci.Update(lset)
ba_vol = [ci.GetElementsOfType(dt,VOL) for dt in [NEG,POS,IF]]
ba_bnd = [ci.GetElementsOfType(dt,BND) for dt in [NEG,POS,IF]]
for dt in [NEG,POS,IF]:
    print(dt, ci.GetElementsOfType(dt,VOL))
for dt in [NEG,POS,IF]:
    print(dt, ci.GetElementsOfType(dt,BND))
# print(ci.GetCutRatios(VOL))
# print(ci.GetCutRatios(BND))

hasneg = BitArray(ci.GetElementsOfType(NEG,VOL))
hasneg |= ci.GetElementsOfType(IF,VOL)
print("hasneg:", hasneg)
haspos = BitArray(ci.GetElementsOfType(POS,VOL))
haspos |= ci.GetElementsOfType(IF,VOL)
print("haspos:", haspos)
has = BitArray(haspos)
has |= hasneg
print("has:", has)

for dt in [NEG,POS,IF]:
    print(dt, ci.GetElementsOfType(dt,VOL))
for dt in [NEG,POS,IF]:
    print(dt, ci.GetElementsOfType(dt,BND))

# gpfacet = GetFacetsWithNeighborTypes(mesh,hasneg,ci.GetElementsOfType(NEG,VOL))
# print(gpfacet)

# gpfacet = GetFacetsWithNeighborTypes(mesh,has,has)
# print(gpfacet)

gpfacet = GetFacetsWithNeighborTypes(mesh,hasneg,ci.GetElementsOfType(IF,VOL))
print(gpfacet)

Draw(lset,mesh,"lset")
