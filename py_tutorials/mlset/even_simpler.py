"""
example of a square described by two level set functions living completely on one element
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
square.AddRectangle([0,0], [1,1], bc=1)
mesh = Mesh (square.GenerateMesh(maxh=2, quad_dominated=False))

lsetp1_a = GridFunction(H1(mesh,order=1))
lsetp1_b = GridFunction(H1(mesh,order=1))

InterpolateToP1( x - 0.5,lsetp1_a)
InterpolateToP1( y - 0.5,lsetp1_b)

Draw (lsetp1_a, mesh, "lset_a")
Draw (lsetp1_b, mesh, "lset_b")

Draw (lsetp1_a * lsetp1_b, mesh, "lset_mult")

##First test: A little square on the right bottom

area = Integrate(levelset_domain = {"levelset": [lsetp1_a,lsetp1_b], "domain_type" : [POS,NEG]},
                 mesh=mesh, cf=1, order=0)

print("area = {:10.8f}".format(area))                     
print("area error = {:4.3e}".format(abs(area-0.25)))     
# inner = { "levelsets" : (lsetp1_a,lsetp1_b,lsetp1_c),
#           "domain_type" : (NEG,NEG,NEG)}

# boundary = { "levelsets" : (lsetp1_a,lsetp1_b,lsetp1_c),
#              "domain_type" : (IF,NEG,NEG) | (NEG,IF,NEG) | (NEG,NEG,IF)}

# outer = { "levelsets" : (lsetp1_a,lsetp1_b,lsetp1_c),
#           "domain_type" : ~(NEG,NEG,NEG)}

# triangle_area = Integrate( levelset_domain=inner, cf=1, mesh=mesh, order=2)
# triangle_boundary_length = Integrate( levelset_domain=boundary, cf=1, mesh=mesh, order=2)
# outer_area = Integrate( levelset_domain=outer, cf=1, mesh=mesh, order=2)

##Second test: A line between (0.5,0.5) and (1,0.5)

length = Integrate(levelset_domain = {"levelset": [lsetp1_a,lsetp1_b], "domain_type" : [POS,IF]},
                 mesh=mesh, cf=1, order=0)

print("length =", length)                     
print("length error =", abs(length-0.5))
