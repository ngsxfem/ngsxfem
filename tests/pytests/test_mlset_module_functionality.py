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

import collections


def CompList(list1, list2):
    return collections.Counter(list1) == collections.Counter(list2)

# -----------------------------------------------------------------------------
# ------------------------------ TEST TRIANGLE --------------------------------
# -----------------------------------------------------------------------------
def test_trainge():
    nr_ls = 3
    triangle = DomainTypeArray([(NEG, NEG, NEG)])

    invert_traingle = ~triangle
    target_list = [(NEG, NEG, POS), (NEG, POS, NEG), (NEG, POS, POS), 
                   (POS, NEG, NEG), (POS, NEG, POS), (POS, POS, NEG), 
                   (POS, POS, POS)]
    assert CompList(invert_traingle.as_list, target_list)

    triangle_bnd = triangle.Boundary()
    target_list = [(IF, NEG, NEG), (NEG, IF, NEG), (NEG, NEG, IF)]
    assert CompList(triangle_bnd.as_list, target_list)
    assert triangle_bnd.codim == 1

    triangle_bnd_bnd = triangle_bnd.Boundary()
    target_list = [(IF, IF, NEG), (NEG, IF, IF), (IF, NEG, IF)]
    assert CompList(triangle_bnd_bnd.as_list, target_list)
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

    target_list = [(POS, NEG), (POS, POS)]
    assert CompList(dta3.as_list, target_list)

    dta3 |= DomainTypeArray([(NEG, POS)])
    target_list = [(POS, NEG), (POS, POS), (NEG, POS)]
    assert CompList(dta3.as_list, target_list)

    dta4 = dta3 & dta1

    target_list = [(POS, NEG)]
    assert CompList(dta4.as_list, target_list)

    dta3 &= dta1 | dta2
    target_list = [(POS, NEG), (POS, POS)]
    assert CompList(dta3.as_list, target_list)

    dta5 = DomainTypeArray([(IF, NEG)]) & DomainTypeArray([(POS, NEG)])
    target_list = [(IF, NEG)]
    assert CompList(dta5.as_list, target_list)

    dta5 &= DomainTypeArray([(POS, IF)])
    target_list = [(IF, IF)]
    assert CompList(dta5.as_list, target_list)

    target_list = [(NEG, POS), (NEG, NEG)]
    assert CompList((~dta3).as_list, target_list)

    dta6 = DomainTypeArray([(IF, NEG, POS)])
    dta6 &= DomainTypeArray([(POS, IF, POS)])
    target_list = [(IF, NEG, IF), (POS, IF, IF), (IF, POS, IF),
                         (NEG, IF, IF), (IF, IF, NEG)]
    assert CompList((~dta6).as_list, target_list)

    dta7 = DomainTypeArray([(NEG, NEG, NEG)])
    target_list = [(NEG, IF, NEG), (NEG, NEG, IF), (IF, NEG, NEG)]
    assert CompList(dta7.Boundary().as_list, target_list)

    dta8 = DomainTypeArray([(POS, IF, NEG)])
    target_list = [(POS, IF, IF), (IF, IF, NEG)]
    assert CompList(dta8.Boundary().as_list, target_list)


# -----------------------------------------------------------------------------
# ---------------------------- TEST ANY EXPANSION -----------------------------
# -----------------------------------------------------------------------------
def test_any_expansion():
    dta9 = DomainTypeArray((NEG, ANY, ANY))
    target_list = [(NEG, NEG, NEG), (NEG, NEG, POS), (NEG, POS, NEG),
                   (NEG, POS, POS)]
    assert CompList(dta9.as_list, target_list)

    dta10 = DomainTypeArray([(NEG, ANY, ANY), (POS, ANY, ANY)])
    dta11 = DomainTypeArray([(ANY, ANY, ANY)])
    assert CompList(dta10.as_list, dta11.as_list)

    dta12 = DomainTypeArray([(IF, ANY, ANY), (POS, IF, ANY)])
    target_list = [(IF, NEG, NEG), (IF, NEG, POS), (IF, POS, NEG), 
                   (IF, POS, POS), (POS, IF, NEG), (POS, IF, POS)]
    assert CompList(dta12.as_list, target_list)


# -----------------------------------------------------------------------------
# --------------------------- RUN TESTS SEPARATELY ----------------------------
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    test_trainge()
    test_operators()
    test_any_expansion()
