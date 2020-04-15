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
mesh = Mesh(geo.GenerateMesh(maxh=1))

# Level sets
V = H1(mesh,order=1)
lsetp1_x_upper, lsetp1_x_lower = GridFunction(V), GridFunction(V)
lsetp1_y_upper, lsetp1_y_lower = GridFunction(V), GridFunction(V)
lsetp1_z_upper, lsetp1_z_lower = GridFunction(V), GridFunction(V)

InterpolateToP1( x - 0.5, lsetp1_x_upper)
InterpolateToP1( x + 0.5, lsetp1_x_lower)
InterpolateToP1( y - 0.5, lsetp1_y_upper)
InterpolateToP1( y + 0.5, lsetp1_y_lower)
InterpolateToP1( z - 0.5, lsetp1_z_upper)
InterpolateToP1( z + 0.5, lsetp1_z_lower)

# Integrate
cube = DomainTypeArray(dtlist=[(NEG, POS, NEG, POS, NEG, POS)])

# Compute volume
volume = Integrate(levelset_domain={"levelset": [lsetp1_x_upper, lsetp1_x_lower, 
                                                 lsetp1_y_upper, lsetp1_y_lower, 
                                                 lsetp1_z_upper, lsetp1_z_lower],
                                    "domain_type": cube.as_list},
                              mesh=mesh, cf=1, order=0)

print("volume = {:10.8f}".format(volume))            
print("volume error = {:4.3e}".format(abs(volume - 1)))
assert abs(volume - 1)


# Compute are of surfaces
for i in range(6):
    domian = [NEG, POS, NEG, POS, NEG, POS]
    domian[i] = IF

    area = Integrate(levelset_domain = {"levelset" : [lsetp1_x_upper, lsetp1_x_lower, 
                                                      lsetp1_y_upper, lsetp1_y_lower, 
                                                      lsetp1_z_upper, lsetp1_z_lower],
                                        "domain_type" : tuple(domian)},
                     mesh=mesh, cf=1, order=0)

    print("area_face{:d} = {:10.8f}".format(i, area))                     
    print("area_face{:d} error = {:4.2e}".format(i, abs(area - 1)))
    assert abs(area-1) < 1e-12

surface_area = Integrate(levelset_domain={"levelset": [lsetp1_x_upper, lsetp1_x_lower, 
                                                 lsetp1_y_upper, lsetp1_y_lower, 
                                                 lsetp1_z_upper, lsetp1_z_lower],
                                    "domain_type": cube.Boundary().as_list},
                              mesh=mesh, cf=1, order=0)

print("surface_area = {:10.8f}".format(surface_area))            
print("surface_area error = {:4.3e}".format(abs(surface_area - 6)))
assert abs(surface_area - 6)


# Compute length of lines
sides = [(IF, POS, IF, POS, NEG, POS), (IF, POS, NEG, IF, NEG, POS),
         (IF, POS, NEG, POS, IF, POS), (IF, POS, NEG, POS, NEG, IF),
         (NEG, IF, IF, POS, NEG, POS), (NEG, IF, NEG, IF, NEG, POS),
         (NEG, IF, NEG, POS, IF, POS), (NEG, IF, NEG, POS, NEG, IF),
         (NEG, POS, IF, POS, IF, POS), (NEG, POS, NEG, IF, IF, POS),
         (NEG, POS, IF, POS, NEG, IF), (NEG, POS, NEG, IF, NEG, IF)]

for i, side in enumerate(sides):
    length = Integrate(levelset_domain = {"levelset": [lsetp1_x_upper, lsetp1_x_lower, 
                                                      lsetp1_y_upper, lsetp1_y_lower, 
                                                      lsetp1_z_upper, lsetp1_z_lower],
                                        "domain_type": side},
                     mesh=mesh, cf=1, order=0)

    
    print("length{:d} = {:10.8f}".format(i, length))                     
    print("length{:d} error = {:4.2e}".format(i, abs(length - 1)))
    assert abs(length - 1) < 1e-12


# Point evaluation
point_domains = {"ppp": (IF, POS, IF, POS, IF, POS),
                 "ppn": (IF, POS, IF, POS, NEG, IF),
                 "pnp": (IF, POS, NEG, IF, IF, POS), 
                 "npp": (NEG, IF, IF, POS, IF, POS),
                 "pnn": (IF, POS, NEG, IF, NEG, IF),
                 "npn": (NEG, IF, IF, POS, NEG, IF),
                 "nnp": (NEG, IF, NEG, IF, IF, POS),
                 "nnn": (NEG, IF, NEG, IF, NEG, IF) }

vals = {"ppp": 1.5, "ppn": 0.5, "pnp": 0.5, "npp": 0.5, "pnn": -0.5, 
        "npn": -0.5, "nnp": -0.5, "nnn": -1.5 }
    
for i, (key, domain) in enumerate(point_domains.items()):

    point = Integrate(levelset_domain = {"levelset": [lsetp1_x_upper, lsetp1_x_lower, 
                                                      lsetp1_y_upper, lsetp1_y_lower, 
                                                      lsetp1_z_upper, lsetp1_z_lower],
                                        "domain_type": domain},
                     mesh=mesh, cf=x+y+z, order=0)

    
    print("point{:d} = {:10.8f}".format(i, point))                     
    print("point{:d} error = {:4.2e}".format(i, abs(point-vals[key])))
    assert abs(point - vals[key]) < 1e-12
