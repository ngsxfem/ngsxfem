"""
	Description:
	------------
	Sharp cube defined by three levelsets.
"""
import math
from netgen.csg import *
from ngsolve import *
from xfem import *


geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-0.8,-0.8,-0.8), Pnt(0.8,0.8,0.8)))
mesh = Mesh(geo.GenerateMesh(maxh=1))


V = H1(mesh,order=1)
lsetp1_x_upper, lsetp1_x_lower, lsetp1_y_upper, lsetp1_y_lower, lsetp1_z_upper, lsetp1_z_lower = GridFunction(V), GridFunction(V), GridFunction(V), GridFunction(V), GridFunction(V), GridFunction(V)

InterpolateToP1( x - 0.5, lsetp1_x_upper)
InterpolateToP1( x + 0.5, lsetp1_x_lower)
InterpolateToP1( y - 0.5, lsetp1_y_upper)
InterpolateToP1( y + 0.5, lsetp1_y_lower)
InterpolateToP1( z - 0.5, lsetp1_z_upper)
InterpolateToP1( z + 0.5, lsetp1_z_lower)

volume = IntegrateMLsetDomain(lsets=[lsetp1_x_upper, lsetp1_x_lower, 
                                     lsetp1_y_upper, lsetp1_y_lower, 
                                     lsetp1_z_upper, lsetp1_z_lower],
                              mesh=mesh,
                              cf=1,
                              order=0,
                              domain_types=[NEG,POS, NEG,POS, NEG, POS])

print("volume = {:10.8f}".format(volume))            
print("volume error = {:4.3e}".format(abs(volume - 1)))

input("")


for i in range(6):
    domian = [NEG,POS, NEG,POS, NEG, POS]
    domian[i] = IF

    area = IntegrateMLsetDomain(lsets=[lsetp1_x_upper, lsetp1_x_lower, 
                                     lsetp1_y_upper, lsetp1_y_lower, 
                                     lsetp1_z_upper, lsetp1_z_lower],
                                mesh=mesh,
                                cf=1,
                                order=0,
                                domain_types=domian)

    print("area{:d} = {:10.8f}".format(i, area))                     
    print("area{:d} error = {:4.2e}".format(i, abs(area-1)))

    input("")   
