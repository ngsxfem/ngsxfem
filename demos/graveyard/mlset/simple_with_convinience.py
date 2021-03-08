"""
Example of a triangle described by three level set functions
"""

from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *

# Background mesh
square = SplineGeometry()
square.AddRectangle([-1,-0.5], [1,1.5], bc=1)
mesh = Mesh (square.GenerateMesh(maxh=1, quad_dominated=False))

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

Draw (IfPos(lsetp1_a,0,1) * IfPos(lsetp1_b,0,1) * IfPos(lsetp1_c,0,1), mesh, "lset_mult")

# Integration over multiple level sets
triangle = DomainTypeArray([(NEG,NEG,NEG)])


area = Integrate(levelset_domain={"levelset": (lsetp1_a,lsetp1_b,lsetp1_c),
                                  "domain_type": triangle.as_list},
                 mesh=mesh, cf=1, order=0)
print("area = {:10.8f}".format(area))                     
print("area error = {:4.3e}".format(abs(area - 0.5)))  
assert abs(area-0.5) < 1e-12

area_invert = Integrate(levelset_domain={"levelset": (lsetp1_a,lsetp1_b,lsetp1_c), 
                                         "domain_type": (~triangle).as_list},
                        mesh=mesh, cf=1, order=0)
print("area_invert = {:10.8f}".format(area_invert))                     
print("area_invert error = {:4.3e}".format(abs(area_invert - 3.5)))  
assert abs(area_invert-3.5) < 1e-12
                   
perimeter = Integrate(levelset_domain={"levelset" : (lsetp1_a,lsetp1_b,lsetp1_c), 
                                       "domain_type": triangle.Boundary().as_list},
                    mesh=mesh, cf=1, order=0)
print("perimeter =", perimeter)                     
print("perimeter error =", abs(perimeter - 1 - sqrt(5)))                     
assert abs(perimeter - 1 - sqrt(5)) < 1e-12


all_points = Integrate(levelset_domain={"levelset": (lsetp1_a,lsetp1_b,lsetp1_c), 
                                        "domain_type": triangle.Boundary().Boundary().as_list},
                        mesh=mesh, cf=y, order=0)
print("all_points =", all_points)                     
print("all_points error =", abs(all_points - 2))
assert abs(all_points - 2) < 1e-12

