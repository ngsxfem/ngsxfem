"""
Test DomainTypeArray functionality
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
import pytest
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *

ngsglobals.msg_level = 1


# -----------------------------------------------------------------------------
# ------------------------------ TEST TRIANGLE --------------------------------
# -----------------------------------------------------------------------------
def test_trianlge():
    nr_ls = 3
    triangle = DomainTypeArray([(NEG, NEG, NEG)])

    invert_traingle = ~triangle
    target = DomainTypeArray([(NEG, NEG, POS), (NEG, POS, NEG), (NEG, POS, POS),
                              (POS, NEG, NEG), (POS, NEG, POS), (POS, POS, NEG), 
                              (POS, POS, POS)])
    assert invert_traingle == target

    triangle_bnd = triangle.Boundary()
    target = DomainTypeArray([(IF, NEG, NEG), (NEG, IF, NEG), (NEG, NEG, IF)])
    assert triangle_bnd == target
    assert triangle_bnd.codim == 1

    target = target | triangle_bnd
    assert target == triangle_bnd

    target = target & triangle_bnd
    assert target == triangle_bnd

    triangle_bnd_bnd = triangle_bnd.Boundary()
    target = DomainTypeArray([(IF, IF, NEG), (NEG, IF, IF), (IF, NEG, IF)])
    assert triangle_bnd_bnd == target
    assert triangle_bnd_bnd.codim == 2

    lines = []
    for i in range(nr_ls):
        line_i = tuple(IF if ii == i else NEG for ii in range(nr_ls))
        lines.append(DomainTypeArray(line_i))
        assert lines[i].codim == 1


# -----------------------------------------------------------------------------
# ------------------------------ TEST OPERATORS -------------------------------
# -----------------------------------------------------------------------------
def test_operators():
    dta1 = DomainTypeArray((POS, NEG))
    dta2 = DomainTypeArray((POS, POS))
    dta3 = dta1 | dta2

    assert (POS, NEG) in dta3
    assert (POS, POS) in dta3
    assert len(dta3) == 2

    target = DomainTypeArray([(POS, NEG), (POS, POS)])
    assert dta3 == target

    dta3 |= DomainTypeArray([(NEG, POS)])
    target = DomainTypeArray([(POS, NEG), (POS, POS), (NEG, POS)])
    assert dta3 == target

    dta4 = dta3 & dta1

    target = DomainTypeArray([(POS, NEG)])
    assert dta4 == target

    dta3 &= dta1 | dta2
    target = DomainTypeArray([(POS, NEG), (POS, POS)])
    assert dta3 == target

    dta5 = DomainTypeArray([(IF, NEG)]) & DomainTypeArray([(POS, NEG)])
    target = DomainTypeArray([(IF, NEG)])
    assert dta5 == target

    dta5 &= DomainTypeArray([(POS, IF)])
    target = DomainTypeArray([(IF, IF)])
    assert dta5 == target

    target = DomainTypeArray([(NEG, POS), (NEG, NEG)])
    assert (~dta3) == target

    dta6 = DomainTypeArray([(IF, NEG, POS)])
    dta6 &= DomainTypeArray([(POS, IF, POS)])
    target = DomainTypeArray([(IF, NEG, IF), (POS, IF, IF), (IF, POS, IF),
                              (NEG, IF, IF), (IF, IF, NEG)])
    assert (~dta6) == target

    dta7 = DomainTypeArray([(NEG, NEG, NEG)])
    target = DomainTypeArray([(NEG, IF, NEG), (NEG, NEG, IF), (IF, NEG, NEG)])
    assert dta7.Boundary() == target

    dta8 = DomainTypeArray([(POS, IF, NEG)])
    target = DomainTypeArray([(POS, IF, IF), (IF, IF, NEG)])
    assert dta8.Boundary() == target


# -----------------------------------------------------------------------------
# ---------------------------- TEST ANY EXPANSION -----------------------------
# -----------------------------------------------------------------------------
def test_any_expansion():
    dta9 = DomainTypeArray((NEG, ANY, ANY))
    target = DomainTypeArray([(NEG, NEG, NEG), (NEG, NEG, POS), 
                            (NEG, POS, NEG), (NEG, POS, POS)])
    assert dta9 == target

    dta10 = DomainTypeArray([(NEG, ANY, ANY), (POS, ANY, ANY)])
    dta11 = DomainTypeArray([(ANY, ANY, ANY)])
    assert dta10 == dta11

    dta12 = DomainTypeArray([(IF, ANY, ANY), (POS, IF, ANY)])
    target = DomainTypeArray([(IF, NEG, NEG), (IF, NEG, POS), (IF, POS, NEG), 
                              (IF, POS, POS), (POS, IF, NEG), (POS, IF, POS)])
    assert dta12 == target


# -----------------------------------------------------------------------------
# -------------------------- TEST TENSOR OPERATORS ----------------------------
# -----------------------------------------------------------------------------
def test_tensor_operators():
    dta1 = DomainTypeArray((POS,))
    dta2 = DomainTypeArray([(POS, NEG), (POS,POS)])
    dta3 = DomainTypeArray((NEG, NEG, IF))

    union = TensorUnion(dta1, dta2)
    target = DomainTypeArray([(POS, ANY, ANY), (ANY, POS, NEG), (ANY, POS, POS)])
    assert union == target

    with pytest.raises(Exception):
        union = TensorUnion(dta1, dta2, dta3)

    intersection1 = TensorIntersection(dta1, dta2)
    target = DomainTypeArray([(POS, POS, NEG), (POS, POS, POS)])
    assert intersection1 == target

    intersection2 = TensorIntersection(dta1, dta2, dta3)
    target = DomainTypeArray([(POS, POS, NEG, NEG, NEG, IF), 
                              (POS, POS,POS, NEG, NEG, IF)])
    assert intersection2 == target


# -----------------------------------------------------------------------------
# --------------------------- RUN TESTS SEPARATELY ----------------------------
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    test_trianlge()
    test_operators()
    test_any_expansion()
    test_tensor_operators()
