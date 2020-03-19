"""
	Description:
	------------
	Sharp cube defined by three levelsets.
"""
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

input("Codim 0")
volume = IntegrateMLsetDomain(lsets=[lsetp1_x_upper, lsetp1_x_lower, 
                                     lsetp1_y_upper, lsetp1_y_lower, 
                                     lsetp1_z_upper, lsetp1_z_lower],
                              mesh=mesh,
                              cf=1,
                              order=0,
                              domain_types=[NEG, POS, NEG, POS, NEG, POS])

print("volume = {:10.8f}".format(volume))            
print("volume error = {:4.3e}".format(abs(volume - 1)))

input("Now codim 1")

for i in range(6):
    domian = [NEG, POS, NEG, POS, NEG, POS]
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

input("Now codim 3")


point_domains = {"ppp": [IF, POS, IF, POS, IF, POS],
                 "ppn": [IF, POS, IF, POS, NEG, IF],
                 "pnp": [IF, POS, NEG, IF, IF, POS], 
                 "npp": [NEG, IF, IF, POS, IF, POS],
                 "pnn": [IF, POS, NEG, IF, NEG, IF],
                 "npn": [NEG, IF, IF, POS, NEG, IF],
                 "nnp": [NEG, IF, NEG, IF, NEG, IF],
                 "nnn": [NEG, IF, NEG, IF, NEG, IF] }

vals = {"ppp": 1.5,
        "ppn": 0.5,
        "pnp": 0.5,
        "npp": 0.5,
        "pnn": -0.5,
        "npn": -0.5,
        "nnp": -0.5,
        "nnn": -1.5
}
    


for key, domain in point_domains.items():

    point = IntegrateMLsetDomain(lsets=[lsetp1_x_upper, lsetp1_x_lower, 
                                     lsetp1_y_upper, lsetp1_y_lower, 
                                     lsetp1_z_upper, lsetp1_z_lower],
                                mesh=mesh,
                                cf=x+y+z,
                                order=0,
                                domain_types=domain)

    print("point{:d} = {:10.8f}".format(i, point))                     
    print("point{:d} error = {:4.2e}".format(i, abs(point-vals[key])))
    input("")