"""
Basic functionality test for multiple level-set functionality
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *

# -----------------------------------------------------------------------------
# --------------------------------- 2D TESTS ----------------------------------
# -----------------------------------------------------------------------------

# ------------------------------ Background Mesh ------------------------------
square = SplineGeometry()
square.AddRectangle([-1, -0.5], [1, 1.5], bc=1)
mesh = Mesh(square.GenerateMesh(maxh=0.5, quad_dominated=False))

# -------------------------------- Level Sets ---------------------------------
level_sets = [y - 1, 2 * x - y, -2 * x - y]
nr_ls = len(level_sets)

level_sets_p1 = tuple(GridFunction(H1(mesh, order=1)) for i in range(nr_ls))
for i, lset_p1 in enumerate(level_sets_p1):
    InterpolateToP1(level_sets[i], lset_p1)

lset_mult = CoefficientFunction(1)
for lset_p1 in level_sets_p1:
    lset_mult *= IfPos(lset_p1, 0, 1)
Draw(lset_mult, mesh, "lset_mult")

# ------------------------------ DomainTypeArray ------------------------------
triangle = DomainTypeArray([(NEG, NEG, NEG)])
Draw(dta_indicator(level_sets_p1, triangle), mesh, "dta_indicator")

# --------------------------------- Cut-Info ----------------------------------
mlci = MultiLevelsetCutInfo(mesh, level_sets_p1)

# --------------------------- Element Marker Tests ----------------------------


def UpdateMarkers(arr, arr_union, arr_intersection=None):
    arr[:] = False
    arr |= arr_union
    if arr_intersection:
        arr &= arr_intersection


def MarkersEqual(arr1, arr2):
    if len(arr1) != len(arr2):
        return False
    for i in range(len(arr1)):
        if arr1[i] != arr2[i]:
            return False
    return True


els_neg, els_hasneg = BitArray(mesh.ne), BitArray(mesh.ne)
els_not_neg, els_if = BitArray(mesh.ne), BitArray(mesh.ne)

UpdateMarkers(els_neg, mlci.GetElementsOfType(triangle.dtlist))
UpdateMarkers(els_hasneg, mlci.GetElementsWithContribution(triangle.dtlist))
UpdateMarkers(els_if, mlci.GetElementsWithContribution(
    triangle.Boundary().dtlist))
UpdateMarkers(
    els_not_neg, mlci.GetElementsWithContribution((~triangle).dtlist))

assert MarkersEqual(~els_neg, els_not_neg)
assert MarkersEqual(els_if, els_hasneg & ~els_neg)

# ----------------------------- Test Integration ------------------------------

# Codim = 0
lset_dom_0 = {"levelset": level_sets_p1[0], "domain_type": NEG}
area0 = Integrate(levelset_domain=lset_dom_0, mesh=mesh, cf=1, order=0)
assert abs(area0 - 3) < 1e-12

lset_dom_tri = {"levelset": level_sets_p1, "domain_type": triangle.dtlist}
area_tri = Integrate(levelset_domain=lset_dom_tri, mesh=mesh, cf=1, order=0)
assert abs(area_tri - 0.5) < 1e-12

lset_dom_tri_i = {"levelset": level_sets_p1, "domain_type": (~triangle).dtlist}
area_tri_i = Integrate(levelset_domain=lset_dom_tri_i,
                       mesh=mesh, cf=1, order=0)
assert abs(area_tri_i - 3.5) < 1e-12

# Codim = 1
lset_dom_side_t = {"levelset": level_sets_p1, "domain_type": (IF, NEG, NEG)}
length_t = Integrate(levelset_domain=lset_dom_side_t, mesh=mesh, cf=1, order=0)
assert abs(length_t - 1) < 1e-12

lset_dom_side_r = {"levelset": level_sets_p1, "domain_type": (NEG, IF, NEG)}
length_r = Integrate(levelset_domain=lset_dom_side_r, mesh=mesh, cf=1, order=0)
assert abs(length_r - sqrt(5 / 4)) < 1e-12

lset_dom_side_l = {"levelset": level_sets_p1, "domain_type": (NEG, NEG, IF)}
length_l = Integrate(levelset_domain=lset_dom_side_l, mesh=mesh, cf=1, order=0)
assert abs(length_l - sqrt(5 / 4)) < 1e-12

lset_dom_per = {"levelset": level_sets_p1,
                "domain_type": triangle.Boundary().dtlist}
length_per = Integrate(levelset_domain=lset_dom_per, mesh=mesh, cf=1, order=0)
assert abs(length_per - 1 - sqrt(5)) < 1e-12

# Codim = 2
lset_dom_pnt_tl = {"levelset": level_sets_p1, "domain_type": (IF, IF, NEG)}
point_val_tl = Integrate(levelset_domain=lset_dom_pnt_tl, mesh=mesh, cf=x + y,
                         order=0)
assert abs(point_val_tl - 1.5) < 1e-12

lset_dom_pnt_tr = {"levelset": level_sets_p1, "domain_type": (IF, NEG, IF)}
point_val_tr = Integrate(levelset_domain=lset_dom_pnt_tr, mesh=mesh, cf=x + y,
                         order=0)
assert abs(point_val_tr - 0.5) < 1e-12

lset_dom_pnt_b = {"levelset": level_sets_p1, "domain_type": (NEG, IF, IF)}
point_val_b = Integrate(levelset_domain=lset_dom_pnt_b, mesh=mesh, cf=x + y,
                        order=0)
assert abs(point_val_b) < 1e-12

lset_dom_cnrs = {"levelset": level_sets_p1,
                 "domain_type": triangle.Boundary().Boundary().dtlist}
point_val_cnrs = Integrate(levelset_domain=lset_dom_cnrs, mesh=mesh, cf=x + y,
                           order=0)
assert abs(point_val_cnrs - 2.0) < 1e-12


del mesh, level_sets, level_sets_p1, mlci, triangle
del els_neg, els_hasneg, els_if, els_not_neg


# -----------------------------------------------------------------------------
# --------------------------------- 3D TESTS ----------------------------------
# -----------------------------------------------------------------------------
from netgen.csg import *
import collections


def CompList(list1, list2):
    return collections.Counter(list1) == collections.Counter(list2)


# ------------------------------ Background Mesh ------------------------------
geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-0.8, -0.8, -0.8), Pnt(0.8, 0.8, 0.8)))
mesh = Mesh(geo.GenerateMesh(maxh=0.5))

# -------------------------------- Level Sets ---------------------------------
level_sets = [x - 0.5, x + 0.5, y - 0.5, y + 0.5, z - 0.5, z + 0.5]
nr_ls = len(level_sets)
level_sets_p1 = tuple(GridFunction(H1(mesh, order=1)) for i in range(nr_ls))
for i, lset_p1 in enumerate(level_sets_p1):
    InterpolateToP1(level_sets[i], lset_p1)

# ------------------------------ DomainTypeArray ------------------------------
cube = DomainTypeArray(dtlist=[(NEG, POS, NEG, POS, NEG, POS)])
sides = []
for i in range(nr_ls):
    domain = [NEG, POS, NEG, POS, NEG, POS]
    domain[i] = IF
    sides.append(DomainTypeArray(tuple(domain)))

target_list = [dtt for side in sides for dtt in side.dtlist]
assert CompList(cube.Boundary().dtlist, target_list)

# --------------------------------- Cut-Info ----------------------------------
mlci = MultiLevelsetCutInfo(mesh, level_sets_p1)

# --------------------------- Element Marker Tests ----------------------------
els_neg, els_hasneg = BitArray(mesh.ne), BitArray(mesh.ne)
els_not_neg, els_if = BitArray(mesh.ne), BitArray(mesh.ne)

UpdateMarkers(els_neg, mlci.GetElementsOfType(cube.dtlist))
UpdateMarkers(els_hasneg, mlci.GetElementsWithContribution(cube.dtlist))
UpdateMarkers(els_if,
              mlci.GetElementsWithContribution(cube.Boundary().dtlist))
UpdateMarkers(els_not_neg,
              mlci.GetElementsWithContribution((~cube).dtlist))

assert MarkersEqual(~els_neg, els_not_neg)
assert MarkersEqual(els_if, els_hasneg & ~els_neg)

els_sides = BitArray(mesh.ne)
els_sides[:] = False
for side in sides:
    els_sides |= mlci.GetElementsWithContribution(side.dtlist)

assert MarkersEqual(els_sides, els_if)

# ----------------------------- Test Integration ------------------------------

# Codim = 0
lset_dom_cube = {"levelset": level_sets_p1, "domain_type": cube.dtlist}
volume_cube = Integrate(levelset_domain=lset_dom_cube,
                        mesh=mesh, cf=1, order=0)
assert abs(volume_cube - 1)

lset_dom_cube_inv = {"levelset": level_sets_p1, "domain_type": (~cube).dtlist}
volume_cube_inv = Integrate(levelset_domain=lset_dom_cube_inv, mesh=mesh, cf=1,
                            order=0)
assert abs(volume_cube_inv - 3.096)

# Codim = 1
for side in sides:
    lset_dom_side = {"levelset": level_sets_p1, "domain_type": side.dtlist}
    area_side = Integrate(levelset_domain=lset_dom_side, mesh=mesh, cf=1,
                          order=0)
    assert abs(area_side - 1) < 1e-12

lset_dom_surface = {"levelset": level_sets_p1,
                    "domain_type": cube.Boundary().dtlist}
area_surface = Integrate(levelset_domain=lset_dom_surface, mesh=mesh, cf=1,
                         order=0)
assert abs(area_surface - 6) < 1e-12

# Codim = 2 is not implemented yet

# Codim = 3

point_domains = {"ppp": (IF, POS, IF, POS, IF, POS),
                 "ppn": (IF, POS, IF, POS, NEG, IF),
                 "pnp": (IF, POS, NEG, IF, IF, POS),
                 "npp": (NEG, IF, IF, POS, IF, POS),
                 "pnn": (IF, POS, NEG, IF, NEG, IF),
                 "npn": (NEG, IF, IF, POS, NEG, IF),
                 "nnp": (NEG, IF, NEG, IF, IF, POS),
                 "nnn": (NEG, IF, NEG, IF, NEG, IF)}

vals = {"ppp": 1.5, "ppn": 0.5, "pnp": 0.5, "npp": 0.5, "pnn": -0.5,
        "npn": -0.5, "nnp": -0.5, "nnn": -1.5}

for i, (key, domain) in enumerate(point_domains.items()):
    lset_pnt = {"levelset": level_sets_p1, "domain_type": domain}
    point = Integrate(levelset_domain=lset_pnt, mesh=mesh, cf=x + y + z,
                      order=0)
    assert abs(point - vals[key]) < 1e-12
