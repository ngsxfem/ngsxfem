"""
example of a triangle described by three level set functions
"""

# the constant pi
from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
# basic geometry features (for the background mesh)
from netgen.geom2d import SplineGeometry

# We generate the background mesh of the domain and use a simplicial triangulation
# To obtain a mesh with quadrilaterals use 'quad_dominated=True'

square = SplineGeometry()
square.AddRectangle([-1,-0.5], [1,1.5], bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.2, quad_dominated=False))

lsetp1_a = GridFunction(H1(mesh,order=1))
lsetp1_b = GridFunction(H1(mesh,order=1))
lsetp1_c = GridFunction(H1(mesh,order=1))

InterpolateToP1(     y-1,lsetp1_a)
InterpolateToP1( 2*x-y  ,lsetp1_b)
InterpolateToP1(-2*x-y  ,lsetp1_c)

Draw (lsetp1_a, mesh, "lset_a")
Draw (lsetp1_b, mesh, "lset_b")
Draw (lsetp1_c, mesh, "lset_c")

Draw (lsetp1_a * lsetp1_b * lsetp1_c, mesh, "lset_mult")

mlci = MultiLevelsetCutInfo(mesh,[lsetp1_a, lsetp1_b, lsetp1_c])
hasneg = mlci.GetElementsOfType([NEG,NEG,NEG])

Draw(BitArrayCF(hasneg),mesh,"hasneg")

hasif = BitArray(mesh.ne)
hasif[:] = False
hasif |= mlci.GetElementsOfType([IF,NEG,NEG])
hasif |= mlci.GetElementsOfType([NEG,IF,NEG])
hasif |= mlci.GetElementsOfType([NEG,NEG,IF])
hasif |= mlci.GetElementsOfType([NEG,IF,IF])
hasif |= mlci.GetElementsOfType([IF,IF,NEG])
hasif |= mlci.GetElementsOfType([IF,NEG,IF])
#later the last five lines should be replaced with:
#hasif = mlci.GetElementsOfType([[NEG,IF,NEG],[IF,NEG,NEG],[NEG,NEG,IF]])

Draw(BitArrayCF(hasif),mesh,"hasif")



