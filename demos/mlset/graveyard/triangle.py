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
square.AddRectangle([-1,-0.5],[1,1.5],bc=1)
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

inner = { "levelset" : (lsetp1_a,lsetp1_b,lsetp1_c),
          "domain_type" : (NEG,NEG,NEG)}

boundary = { "levelset" : (lsetp1_a,lsetp1_b,lsetp1_c),
             "domain_type" : [(IF,NEG,NEG), (NEG,IF,NEG), (NEG,NEG,IF)]}

# outer = { "levelsets" : (lsetp1_a,lsetp1_b,lsetp1_c),
#           "domain_type" : ~(NEG,NEG,NEG)}

triangle_area = Integrate( levelset_domain=inner, cf=1, mesh=mesh, order=2)
triangle_boundary_length = Integrate( levelset_domain=boundary, cf=1, mesh=mesh, order=2)
# outer_area = Integrate( levelset_domain=outer, cf=1, mesh=mesh, order=2)
