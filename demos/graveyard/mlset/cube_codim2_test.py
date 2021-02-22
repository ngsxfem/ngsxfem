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

# -----------------------------------------------------------------------------
# First Test

# Level sets
level_sets = (x - 0.5, x + 0.5, x - y, z - 0)
nr_ls = len(level_sets)
level_sets_p1 = tuple(GridFunction(H1(mesh,order=1)) for i in range(nr_ls))

for i, lset_p1 in enumerate(level_sets_p1):
    InterpolateToP1(level_sets[i], lset_p1)

# Define Domain
line = DomainTypeArray(dtlist=[(NEG, POS, IF, IF)])
Draw(line.IndicatorSmoothed(level_sets_p1), mesh, "line", sd=6)

# Compute length
length = Integrate(levelset_domain={"levelset": level_sets_p1,
                                    "domain_type": line},
                   mesh=mesh, cf=1, order=0)

print("length = {:10.8f}".format(length))
print("length error = {:4.3e}".format(abs(length - sqrt(2))))
assert abs(length - sqrt(2)) < 1e-12

del level_sets, level_sets_p1, nr_ls, line

# -----------------------------------------------------------------------------
# Second Test: cross

# Level sets
level_sets = (x - 0.5, - x - 0.5, z - y, x - z, x + z)
nr_ls = len(level_sets)
level_sets_p1 = tuple(GridFunction(H1(mesh,order=1)) for i in range(nr_ls))

for i, lset_p1 in enumerate(level_sets_p1):
    InterpolateToP1(level_sets[i], lset_p1)

# Define Domain
cross = DomainTypeArray([(NEG, NEG, IF, IF, ANY), (NEG, NEG, IF, ANY, IF)])
Draw(cross.IndicatorSmoothed(level_sets_p1), mesh, "cross", sd=6)

# Compute length
length = Integrate(levelset_domain={"levelset": level_sets_p1,
                                    "domain_type": cross},
                   mesh=mesh, cf=1, order=0)

print("length = {:10.8f}".format(length))
print("length error = {:4.3e}".format(abs(length - 2*sqrt(3))))
assert abs(length - 2*sqrt(3)) < 1e-12

del level_sets, level_sets_p1, nr_ls, cross

# -----------------------------------------------------------------------------
# Third test: Parallel lines close to each other

level_sets = (x - 0.5, - x - 0.5, - z - 0.5, z + y, x + y, z + y + 0.1)
nr_ls = len(level_sets)
level_sets_p1 = tuple(GridFunction(H1(mesh,order=1)) for i in range(nr_ls))

for i, lset_p1 in enumerate(level_sets_p1):
    InterpolateToP1(level_sets[i], lset_p1)

lines = DomainTypeArray([(NEG, NEG, NEG, IF, IF, ANY), 
                         (NEG, NEG, NEG, ANY, IF, IF)])

length = Integrate(levelset_domain={"levelset": level_sets_p1,
                                    "domain_type": lines},
                   mesh=mesh, cf=1, order=0)

print("length = {:10.8f}".format(length))
print("length error = {:4.3e}".format(abs(length - 1.9*sqrt(3))))
assert abs(length - 1.9*sqrt(3)) < 1e-12