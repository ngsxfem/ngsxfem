"""
	Description:
	------------
	Sharp cube defined by three levelsets.
"""

# Import libraries
from netgen.csg import *
from ngsolve import *
from xfem import *
from xfem.mlset import *

# Mesh
geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-0.8,-0.8,-0.8), Pnt(0.8,0.8,0.8)))
mesh = Mesh(geo.GenerateMesh(maxh=0.75))

# Level set
V = H1(mesh,order=1)
lsetp1_1, lsetp1_2, lsetp1_3 = GridFunction(V), GridFunction(V), GridFunction(V)

def cf_abs(x):
	return IfPos(x,x,-x)
InterpolateToP1( cf_abs(x) - 0.5, lsetp1_1)
InterpolateToP1( cf_abs(y) - 0.5, lsetp1_2)
InterpolateToP1( cf_abs(z) - 0.5, lsetp1_3)


Draw(lsetp1_1,mesh,"lsetp1_1")
Draw(lsetp1_2,mesh,"lsetp1_2")
Draw(lsetp1_3,mesh,"lsetp1_3")

# Integrate
cube = DomainTypeArray(dtlist=[(NEG,NEG,NEG)])

volume = Integrate(levelset_domain={"levelset": [lsetp1_1, lsetp1_2, lsetp1_3], 
                                    "domain_type": cube.as_list},
                   mesh=mesh, cf=1, order=0)


print("volume = {:10.8f}".format(volume))            
print("volume error = {:4.3e}".format(abs(volume - 1)))  
assert abs(volume - 1) < 1e-12

surface_area = Integrate(levelset_domain={"levelset": [lsetp1_1, lsetp1_2, lsetp1_3], 
                                          "domain_type": cube.Boundary().as_list},
                         mesh=mesh, cf=1, order=0)

print("surface_area = {:10.8f}".format(surface_area))            
print("surface_area error = {:4.3e}".format(abs(surface_area - 6)))  
assert abs(surface_area - 6) < 1e-12