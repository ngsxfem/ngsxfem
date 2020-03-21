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
mesh = Mesh (square.GenerateMesh(maxh=0.3, quad_dominated=False))

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

area0 = Integrate(levelset_domain={"levelset" : lsetp1_a, "domain_type": NEG},
                  mesh=mesh, cf=1, order=0)
print("area0 = {:10.8f}".format(area0))                     
print("area error = {:4.3e}".format(abs(area0-3)))

input("")

area = Integrate(levelset_domain={"levelset" : [lsetp1_a,lsetp1_b,lsetp1_c], "domain_type": [NEG,NEG,NEG]},
                         mesh=mesh, cf=1, order=0)
print("area = {:10.8f}".format(area))                     
print("area error = {:4.3e}".format(abs(area-0.5)))  
                   
input("")

area = Integrate(levelset_domain={"levelset" : [lsetp1_a,lsetp1_b,lsetp1_c], "domain_type": [[NEG,NEG,NEG]]},
                         mesh=mesh, cf=1, order=0)
print("[duplicate] area = {:10.8f}".format(area))                     
print("[duplicate] area error = {:4.3e}".format(abs(area-0.5)))  
                   
input("")

length1 = Integrate(levelset_domain={"levelset" : [lsetp1_a,lsetp1_b,lsetp1_c], "domain_type": [IF,NEG,NEG]},
                         mesh=mesh, cf=1, order=0)
print("length1 =", length1)                     
print("length1 error =", abs(length1-1))                     

input("")

length2 = Integrate(levelset_domain={"levelset" : [lsetp1_a,lsetp1_b,lsetp1_c], "domain_type": [NEG,IF,NEG]},
                         mesh=mesh, cf=1, order=0)
print("length2 =", length2)                     
print("length2 error =", abs(length2-sqrt(5/4)))                     

input("")

length3 = Integrate(levelset_domain={"levelset" : [lsetp1_a,lsetp1_b,lsetp1_c], "domain_type": [NEG,NEG,IF]},
                         mesh=mesh, cf=1, order=0)
print("length3 =", length3)                     
print("length3 error =", abs(length3-sqrt(5/4)))                     

input("")
point_val_y = Integrate(levelset_domain={"levelset" : [lsetp1_a,lsetp1_b,lsetp1_c], "domain_type": [NEG,IF,IF]},
                        mesh=mesh, cf=y, order=0)
print("point_val_y =", point_val_y)                     
print("point_val_y error =", abs(point_val_y-0))

input("")

point_val_yp = Integrate(levelset_domain={"levelset" : [lsetp1_a,lsetp1_b,lsetp1_c], "domain_type": [NEG,IF,IF]},
                         mesh=mesh, cf=1-y, order=0)
print("point_val_yp =", point_val_yp)                     
print("point_val_yp error =", abs(point_val_yp-1))

input("")

# inner = { "levelsets" : (lsetp1_a,lsetp1_b,lsetp1_c),
#           "domain_type" : (NEG,NEG,NEG)}

# boundary = { "levelsets" : (lsetp1_a,lsetp1_b,lsetp1_c),
#              "domain_type" : (IF,NEG,NEG) | (NEG,IF,NEG) | (NEG,NEG,IF)}

# outer = { "levelsets" : (lsetp1_a,lsetp1_b,lsetp1_c),
#           "domain_type" : ~(NEG,NEG,NEG)}

# triangle_area = Integrate( levelset_domain=inner, cf=1, mesh=mesh, order=2)
# triangle_boundary_length = Integrate( levelset_domain=boundary, cf=1, mesh=mesh, order=2)
# outer_area = Integrate( levelset_domain=outer, cf=1, mesh=mesh, order=2)
