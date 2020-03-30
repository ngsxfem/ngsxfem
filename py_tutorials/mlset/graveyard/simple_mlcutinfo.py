"""
example of a triangle described by three level set functions
"""

from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *


# Background mesh
square = SplineGeometry()
square.AddRectangle([-1, -0.5], [1, 1.5], bc=1)
mesh = Mesh(square.GenerateMesh(maxh=0.2, quad_dominated=False))

# Level sets
lsetp1_a = GridFunction(H1(mesh, order=1))
lsetp1_b = GridFunction(H1(mesh, order=1))
lsetp1_c = GridFunction(H1(mesh, order=1))

InterpolateToP1(y - 1, lsetp1_a)
InterpolateToP1(2 * x - y, lsetp1_b)
InterpolateToP1(-2 * x - y, lsetp1_c)

Draw(lsetp1_a, mesh, "lset_a")
Draw(lsetp1_b, mesh, "lset_b")
Draw(lsetp1_c, mesh, "lset_c")

# Shape of interest
triangle = DomainTypeArray([(NEG, NEG, NEG)])
Draw(triangle.Indicator((lsetp1_a, lsetp1_b, lsetp1_c)), mesh, "triangle_dom")
Draw(triangle.Boundary().IndicatorSomoothed((lsetp1_a, lsetp1_b, lsetp1_c)), 
                                            mesh, "triangle_bnd")

# Show MultiLevelsetCutInfo functionality
mlci = MultiLevelsetCutInfo(mesh, (lsetp1_a, lsetp1_b, lsetp1_c))

els_neg = mlci.GetElementsOfType(triangle)
els_neg2 = mlci.GetElementsOfType(triangle.as_list)
els_not_neg = mlci.GetElementsOfType(~triangle)
els_if = mlci.GetElementsWithContribution(triangle.Boundary())
els_hasneg = mlci.GetElementsWithContribution(triangle)
els_hasneg2 = mlci.GetElementsWithContribution(triangle.as_list)

els_if2 = BitArray(mesh.ne)
els_if2[:] = False
els_if2 |= els_hasneg & ~els_neg

# Draw BitArrays
Draw(BitArrayCF(els_neg), mesh, "els_neg")
Draw(BitArrayCF(els_neg2), mesh, "els_neg2")
Draw(BitArrayCF(els_hasneg), mesh, "els_hasneg")
Draw(BitArrayCF(els_hasneg2), mesh, "els_hasneg2")
Draw(BitArrayCF(els_if), mesh, "els_if")
Draw(BitArrayCF(els_if2), mesh, "els_if2")
