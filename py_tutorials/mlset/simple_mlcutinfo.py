"""
example of a triangle described by three level set functions
"""

from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *


# Background mesh
square = SplineGeometry()
square.AddRectangle([-1,-0.5], [1,1.5], bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.1, quad_dominated=False))

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

# Shape of interest
triangle = DomainTypeArray([(NEG,NEG,NEG)])

# Show MultiLevelsetCutInfo functionality
mlci = MultiLevelsetCutInfo(mesh,(lsetp1_a, lsetp1_b, lsetp1_c))

els_neg = mlci.GetElementsOfType(triangle.dtlist)
els_not_neg = mlci.GetElementsOfType((~triangle).dtlist)
els_outer = mlci.GetElementsOfType((~triangle).dtlist + (~(triangle.Boundary())).dtlist)
els_if1 = mlci.GetElementsOfType(triangle.Boundary(element_marking=True)) 
els_if2 = mlci.GetElementsOfType(triangle.Boundary(element_marking=False)) 

els_hasneg = BitArray(mesh.ne)
els_hasneg.Clear()
els_hasneg |= els_neg | els_if1

# Draw BitArrays
Draw(BitArrayCF(els_neg), mesh, "els_neg")
Draw(BitArrayCF(els_if1), mesh, "els_if1")
Draw(BitArrayCF(els_if2), mesh, "els_if2")
Draw(BitArrayCF(els_not_neg), mesh, "els_not_neg")
Draw(BitArrayCF(els_outer), mesh, "els_outer")
Draw(BitArrayCF(els_hasneg), mesh, "els_hasneg")


