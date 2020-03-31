"""
    Description:
    -----------
    Zalesak Disk described my multiple level sets. Inspired by [1].

    References:
    -----------
    [1] D. P. Starinshak, S. Karni and P. L. Roe: A New Level-Set
        Model for the Representation of Non-Smooth Geometries. In J Sci 
        Comput (2014) 61(3):649-672. DOI:10.1007/s10915-014-9842-0

"""

from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *
from math import pi, asin, sqrt

SetNumThreads(4)

geo = SplineGeometry()
geo.AddRectangle((-1.1, -1.1), (1.1, 1.1), bc=1)

with TaskManager():
    mesh = Mesh(geo.GenerateMesh(maxh=0.01))

level_sets = [x * x + y * y - 1, -x - 1 / 3, x - 1 / 3, y - 0.5]
nr_ls = len(level_sets)
level_sets_p1 = tuple(GridFunction(H1(mesh, order=1)) for i in range(nr_ls))

for i, lsetp1 in enumerate(level_sets_p1):
    InterpolateToP1(level_sets[i], lsetp1)
    Draw(lsetp1, mesh, "lset_p1_" + str(i))

SetVisualization(min=0, max=0)
Redraw()

area_zd = pi - 2 / 3 * (0.5 + sqrt(2) * 2 / 3)
area_zd += - (2 * asin(1 / 3) - sin(2 * asin(1 / 3))) / 2
print("Area   = {:12.10f}".format(area_zd))


z_disc1 = DomainTypeArray([(NEG, NEG, NEG, POS), (NEG, POS, NEG, POS),
                           (NEG, POS, NEG, NEG), (NEG, NEG, POS, NEG),
                           (NEG, NEG, POS, POS)])
lset_zdisc1 = {"levelset": level_sets_p1, "domain_type": z_disc1}

area1 = Integrate(levelset_domain=lset_zdisc1, cf=1, mesh=mesh, order=0)
print("Area 1 = {:12.10f}".format(area1))
print("Error: {:4.2e}".format(abs(area1 - area_zd)))


# Alternative 1
z_disc2 = DomainTypeArray([(NEG, ANY, ANY, POS), (NEG, POS, ANY, ANY),
                           (NEG, ANY, POS, ANY)])
lset_zdisc2 = {"levelset": level_sets_p1, "domain_type": z_disc2}

area2 = Integrate(levelset_domain=lset_zdisc2, cf=1, mesh=mesh, order=0)
print("Area 2 = {:12.10f}".format(area2))
print("Error: {:4.2e}".format(abs(area2 - area_zd)))
assert abs(area1 - area2) < 1e-12

# # Alternative 2
part1 = DomainTypeArray((NEG,ANY,ANY,ANY))
part2 = DomainTypeArray((NEG,NEG,NEG,NEG))
z_disc3 = part1 & ~part2

lset_zdisc3 = {"levelset": level_sets_p1, "domain_type": z_disc3}

area3 = Integrate(levelset_domain=lset_zdisc3, cf=1, mesh=mesh, order=0)
print("Area 3 = {:12.10f}".format(area3))
print("Error: {:4.2e}".format(abs(area3 - area_zd)))
assert abs(area1 - area3) < 1e-12
