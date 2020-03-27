"""
Basic functionality test for multiple level-set functionality
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
import pytest
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *

ngsglobals.msg_level = 1


# ----------------------------- UTILITY FUNTIONS ------------------------------
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


# -----------------------------------------------------------------------------
# --------------------------------- 2D TESTS ----------------------------------
# -----------------------------------------------------------------------------


def test_2d_mlci_and_lo_integration():
    # ---------------------------- Background Mesh ----------------------------
    square = SplineGeometry()
    square.AddRectangle([-1, -0.5], [1, 1.5], bc=1)
    mesh = Mesh(square.GenerateMesh(maxh=0.5, quad_dominated=False))

    # ------------------------------ Level Sets -------------------------------
    level_sets = [y - 1, 2 * x - y, -2 * x - y]
    nr_ls = len(level_sets)

    level_sets_p1 = tuple(GridFunction(H1(mesh, order=1)) for i in range(nr_ls))
    for i, lset_p1 in enumerate(level_sets_p1):
        InterpolateToP1(level_sets[i], lset_p1)

    lset_mult = CoefficientFunction(1)
    for lset_p1 in level_sets_p1:
        lset_mult *= IfPos(lset_p1, 0, 1)
    Draw(lset_mult, mesh, "lset_mult")

    # ---------------------------- DomainTypeArray ----------------------------
    triangle = DomainTypeArray((NEG, NEG, NEG))
    Draw(dta_indicator(level_sets_p1, triangle), mesh, "dta_indicator")

    # ------------------------------- Cut-Info --------------------------------
    mlci = MultiLevelsetCutInfo(mesh, level_sets_p1)

    # ------------------------- Element Marker Tests --------------------------
    els_neg, els_hasneg = BitArray(mesh.ne), BitArray(mesh.ne)
    els_neg2, els_hasneg2 = BitArray(mesh.ne), BitArray(mesh.ne)
    els_not_neg, els_if = BitArray(mesh.ne), BitArray(mesh.ne)

    UpdateMarkers(els_neg, mlci.GetElementsOfType(triangle))
    UpdateMarkers(els_hasneg, mlci.GetElementsWithContribution(triangle))
    UpdateMarkers(els_neg2, mlci.GetElementsOfType(triangle.as_list))
    UpdateMarkers(els_hasneg2,
                  mlci.GetElementsWithContribution(triangle.as_list))
    UpdateMarkers(els_if, mlci.GetElementsWithContribution(triangle.Boundary()))
    UpdateMarkers(els_not_neg, mlci.GetElementsWithContribution(~triangle))

    assert MarkersEqual(els_neg, els_neg2)
    assert MarkersEqual(els_hasneg, els_hasneg2)
    assert MarkersEqual(~els_neg, els_not_neg)
    assert MarkersEqual(els_if, els_hasneg & ~els_neg)

    # --------------------------- Test Integration ----------------------------

    # ---------------------- Test Low-Order Integration -----------------------
    # Codim = 0
    lset_dom_0 = {"levelset": level_sets_p1[0], "domain_type": NEG}
    area0 = Integrate(levelset_domain=lset_dom_0, mesh=mesh, cf=1, order=0)
    assert abs(area0 - 3) < 1e-12

    lset_dom_tri = {"levelset": level_sets_p1, "domain_type": triangle}
    area_tri = Integrate(levelset_domain=lset_dom_tri, mesh=mesh, cf=1, order=0)
    assert abs(area_tri - 0.5) < 1e-12

    lset_dom_tri_i = {"levelset": level_sets_p1, "domain_type": (~triangle)}
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
                    "domain_type": triangle.Boundary()}
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
                     "domain_type": triangle.Boundary().Boundary()}
    point_val_cnrs = Integrate(levelset_domain=lset_dom_cnrs, mesh=mesh, cf=x + y,
                               order=0)
    assert abs(point_val_cnrs - 2.0) < 1e-12


    del mesh, level_sets, level_sets_p1, mlci, triangle
    del els_neg, els_hasneg, els_if, els_not_neg


def test_2d_ho_integration():
    # -------------------- Test Higher-Order Integration ----------------------
    level_sets = [-y, x - 1, y - 1, -x]
    nr_ls = len(level_sets)

    # ---------------------------- Background Mesh ----------------------------
    geo = SplineGeometry()
    geo.AddRectangle((-0.2, -0.2), (1.2, 1.2),
                     bcs=("bottom", "right", "top", "left"))
    mesh = Mesh(geo.GenerateMesh(maxh=0.4))

    # ------------------------------- LEVELSET --------------------------------
    level_sets_p1 = tuple(GridFunction(H1(mesh, order=1)) for i in range(nr_ls))
    for i, lsetp1 in enumerate(level_sets_p1):
        InterpolateToP1(level_sets[i], lsetp1)

    square = DomainTypeArray([(NEG, NEG, NEG, NEG)])
    lset_dom_square = {"levelset": level_sets_p1, "domain_type": square}

    # ------------------------------- Integrate -------------------------------
    result = Integrate(levelset_domain=lset_dom_square, mesh=mesh, 
                       cf=x * (1 - x) * y * (1 - y), order=4)
    assert abs(result - 1/36) < 1e-12

    del mesh, level_sets, level_sets_p1, square 


# -----------------------------------------------------------------------------
# --------------------------------- 3D TESTS ----------------------------------
# -----------------------------------------------------------------------------
from netgen.csg import *
import collections


def CompList(list1, list2):
    return collections.Counter(list1) == collections.Counter(list2)

def test_3d_mlset():
    # ---------------------------- Background Mesh ----------------------------
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-0.8, -0.8, -0.8), Pnt(0.8, 0.8, 0.8)))
    mesh = Mesh(geo.GenerateMesh(maxh=0.5))

    # ------------------------------ Level Sets -------------------------------
    level_sets = [x - 0.5, x + 0.5, y - 0.5, y + 0.5, z - 0.5, z + 0.5]
    nr_ls = len(level_sets)
    level_sets_p1 = tuple(GridFunction(H1(mesh, order=1)) for i in range(nr_ls))
    for i, lset_p1 in enumerate(level_sets_p1):
        InterpolateToP1(level_sets[i], lset_p1)

    # ---------------------------- DomainTypeArray ----------------------------
    cube = DomainTypeArray(dtlist=[(NEG, POS, NEG, POS, NEG, POS)])
    sides = []
    for i in range(nr_ls):
        domain = [NEG, POS, NEG, POS, NEG, POS]
        domain[i] = IF
        sides.append(DomainTypeArray(tuple(domain)))

    target_list = [dtt for side in sides for dtt in side.as_list]
    assert CompList(cube.Boundary().as_list, target_list)

    # ------------------------------- Cut-Info --------------------------------
    mlci = MultiLevelsetCutInfo(mesh, level_sets_p1)

    # ------------------------- Element Marker Tests --------------------------
    els_neg, els_hasneg = BitArray(mesh.ne), BitArray(mesh.ne)
    els_not_neg, els_if = BitArray(mesh.ne), BitArray(mesh.ne)

    UpdateMarkers(els_neg, mlci.GetElementsOfType(cube))
    UpdateMarkers(els_hasneg, mlci.GetElementsWithContribution(cube))
    UpdateMarkers(els_if, mlci.GetElementsWithContribution(cube.Boundary()))
    UpdateMarkers(els_not_neg, mlci.GetElementsWithContribution(~cube))

    assert MarkersEqual(~els_neg, els_not_neg)
    assert MarkersEqual(els_if, els_hasneg & ~els_neg)

    els_sides = BitArray(mesh.ne)
    els_sides[:] = False
    for side in sides:
        els_sides |= mlci.GetElementsWithContribution(side)

    assert MarkersEqual(els_sides, els_if)

    # --------------------------- Test Integration ----------------------------

    # Codim = 0
    lset_dom_cube = {"levelset": level_sets_p1, "domain_type": cube}
    volume_cube = Integrate(levelset_domain=lset_dom_cube,
                            mesh=mesh, cf=1, order=0)
    assert abs(volume_cube - 1)

    lset_dom_cube_inv = {"levelset": level_sets_p1, "domain_type": ~cube}
    volume_cube_inv = Integrate(levelset_domain=lset_dom_cube_inv, mesh=mesh, cf=1,
                                order=0)
    assert abs(volume_cube_inv - 3.096)

    # Codim = 1
    for side in sides:
        lset_dom_side = {"levelset": level_sets_p1, "domain_type": side}
        area_side = Integrate(levelset_domain=lset_dom_side, mesh=mesh, cf=1,
                              order=0)
        assert abs(area_side - 1) < 1e-12

    lset_dom_surface = {"levelset": level_sets_p1,
                        "domain_type": cube.Boundary()}
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

# -----------------------------------------------------------------------------
# --------------------------- RUN TESTS SEPARATELY ----------------------------
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    test_2d_mlci_and_lo_integration()
    test_2d_ho_integration()
    test_3d_mlset()