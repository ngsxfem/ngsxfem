"""
example of a triangle described by three level set functions
"""

from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *


# Background meshsquare = SplineGeometry()
square.AddRectangle([-1,-0.5], [1,1.5], bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.2, quad_dominated=False))

# Level sets
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

# Show MultiLevelsetCutInfo functionality
mlci = MultiLevelsetCutInfo(mesh,(lsetp1_a, lsetp1_b, lsetp1_c))

hasneg = mlci.GetElementsOfType((NEG,NEG,NEG))

Draw(BitArrayCF(hasneg),mesh,"hasneg")

hasif = mlci.GetElementsOfType([(IF,NEG,NEG),(NEG,IF,NEG),(NEG,NEG,IF),(IF,IF,NEG),(IF,NEG,IF),(NEG,IF,IF)])

Draw(BitArrayCF(hasif),mesh,"hasif")



