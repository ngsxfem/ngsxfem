"""
	Description:
	------------
	Sharp cube defined by six level sets, one for each face.
"""

# Import libraries
from netgen.csg import *
from ngsolve import *
from xfem import *
from xfem.mlset import *

# Mesh
geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-0.8,-0.8,-0.8), Pnt(0.8,0.8,0.8)))
mesh = Mesh(geo.GenerateMesh(maxh=0.3))

# Level sets
V = H1(mesh,order=1)
lsetp1_x_upper, lsetp1_x_lower = GridFunction(V), GridFunction(V)
lsetp1_y = GridFunction(V)
lsetp1_z = GridFunction(V)

InterpolateToP1( x - 0.5, lsetp1_x_upper)
InterpolateToP1( x + 0.5, lsetp1_x_lower)
InterpolateToP1( x - y , lsetp1_y)
InterpolateToP1( z - 0., lsetp1_z)

# Integrate
cube = DomainTypeArray(dtlist=[(NEG, POS, IF, IF)])

# Compute length
length = Integrate(levelset_domain={"levelset": [lsetp1_x_upper, lsetp1_x_lower, 
                                                 lsetp1_y, lsetp1_z],
                                    "domain_type": cube.as_list},
                              mesh=mesh, cf=1, order=0)

print("length = {:10.8f}".format(length))
print("length error = {:4.3e}".format(abs(length - sqrt(2))))
assert abs(length - sqrt(2)) < 1e-12
