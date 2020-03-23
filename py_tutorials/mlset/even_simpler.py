"""
example of a square described by two level set functions living completely on one element
"""

# import netgen and xfem stuff
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *


# Generate mesh
square = SplineGeometry()
square.AddRectangle([0,0], [1,1], bc=1)
mesh = Mesh (square.GenerateMesh(maxh=2, quad_dominated=False))

# Levelsets
lsetp1_a = GridFunction(H1(mesh,order=1))
lsetp1_b = GridFunction(H1(mesh,order=1))

InterpolateToP1( x - 0.5,lsetp1_a)
InterpolateToP1( y - 0.5,lsetp1_b)

Draw (lsetp1_a, mesh, "lset_a")
Draw (lsetp1_b, mesh, "lset_b")
Draw (lsetp1_a * lsetp1_b, mesh, "lset_mult")


# First test: A little square on the right bottom
area = Integrate(levelset_domain={"levelset": [lsetp1_a,lsetp1_b], 
                                  "domain_type": (POS,NEG)},
                 mesh=mesh, cf=1, order=0)

print("area = {:10.8f}".format(area))                     
print("area error = {:4.3e}".format(abs(area-0.25)))     


length1 = Integrate(levelset_domain={"levelset": [lsetp1_a, lsetp1_b], 
                                     "domain_type": (POS,IF)},
                    mesh=mesh, cf=1, order=0)

print("length1 =", length1)                     
print("length1 error =", abs(length1-0.5))
assert abs(length1-0.5) < 1e-12

length2 = Integrate(levelset_domain={"levelset": [lsetp1_a, lsetp1_b], 
                                     "domain_type": [(POS, IF), (IF, NEG)]},
                    mesh=mesh, cf=1, order=0)

print("length2 =", length2)                     
print("length2 error =", abs(length2-1))
assert abs(length2-1) < 1e-12