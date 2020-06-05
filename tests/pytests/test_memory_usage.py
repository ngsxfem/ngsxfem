"""
Test for memory usage/clean-up in integrate
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
import pytest
from netgen.geom2d import unit_square
from ngsolve import *
from xfem import *

import psutil
import os

# ----------------------------- UTILITY FUNTIONS ------------------------------

def memory_usage_psutil():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2 ** 20)
    return mem


# -----------------------------------------------------------------------------
# ------------------------------ SINGLE LEVELSET ------------------------------
# -----------------------------------------------------------------------------
def test_singele_lset_memory():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    
    lset = (x-0.5)**2 + (y - 0.5)**2 - 0.4**2
    lset_p1 = GridFunction(H1(mesh, order=1))
    InterpolateToP1(lset, lset_p1)

    lset_dom = {"levelset": lset_p1, "domain_type": NEG, "subdivions": 0}

    m1 = memory_usage_psutil()

    for it in range(int(1e4)):
        Integrate(lset_dom, cf=1, mesh=mesh, order=0)

    m2 = memory_usage_psutil()
    assert abs(m2 - m1) < 0.4

# -----------------------------------------------------------------------------
# ---------------------------- MULTIPLE LEVELSETS -----------------------------
# -----------------------------------------------------------------------------
def test_multiple_lset_memory():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    
    lsets = [-y + 0.1, x - 0.9, y - 0.9, -x + 0.1]
    lsets_p1 = tuple(GridFunction(H1(mesh, order=1)) for i in range(len(lsets)))
    for lset, lsetp1 in zip(lsets, lsets_p1):
        InterpolateToP1(lset, lsetp1)

    lset_dom = {"levelset": lsets_p1, "domain_type": (NEG, NEG, NEG, NEG)}

    m1 = memory_usage_psutil()

    for it in range(int(1e4)):
        Integrate(lset_dom, cf=1, mesh=mesh, order=0)

    m2 = memory_usage_psutil()
    assert abs(m2 - m1) < 0.4


# -----------------------------------------------------------------------------
# -------------------------------- SPACE-TIME ---------------------------------
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    test_singele_lset_memory()
    test_multiple_lset_memory()