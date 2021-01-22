"""
Test for memory usage/clean-up in integrate
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
import pytest

from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh
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
@pytest.mark.parametrize("quad", [True, False])
def test_singele_lset_memory(quad):
    mesh = MakeStructured2DMesh(quads=quad, nx=2**3, ny=2**3)
    
    lset = (x-0.5)**2 + (y - 0.5)**2 - 0.4**2
    lset_p1 = GridFunction(H1(mesh, order=1))
    InterpolateToP1(lset, lset_p1)

    lset_dom = {"levelset": lset_p1, "domain_type": NEG, "subdivions": 0}

    mem = []
    mem.append(memory_usage_psutil())

    for it in range(2):
        for it in range(int(5e3)):
            Integrate(lset_dom, cf=1, mesh=mesh, order=0)

        mem.append(memory_usage_psutil())

    assert abs(mem[1] - mem[2]) <= abs(max(mem) - min(mem))
    assert abs(mem[1] - mem[2]) < 0.01 * max(mem)


# -----------------------------------------------------------------------------
# ---------------------------- MULTIPLE LEVELSETS -----------------------------
# -----------------------------------------------------------------------------
def test_multiple_lset_memory():
    mesh = MakeStructured2DMesh(quads=False, nx=2**3, ny=2**3)
    
    lsets = [-y + 0.1, x - 0.9, y - 0.9, -x + 0.1]
    lsets_p1 = tuple(GridFunction(H1(mesh, order=1)) for i in range(len(lsets)))
    for lset, lsetp1 in zip(lsets, lsets_p1):
        InterpolateToP1(lset, lsetp1)

    lset_dom = {"levelset": lsets_p1, "domain_type": (NEG, NEG, NEG, NEG)}

    mem = []
    mem.append(memory_usage_psutil())

    for it in range(2):
        for it in range(int(5e3)):
            Integrate(lset_dom, cf=1, mesh=mesh, order=0)
        mem.append(memory_usage_psutil())

    assert abs(mem[1] - mem[2]) <= abs(max(mem) - min(mem))
    assert abs(mem[1] - mem[2]) < 0.01 * max(mem)


# -----------------------------------------------------------------------------
# -------------------------------- SPACE-TIME ---------------------------------
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    test_singele_lset_memory(True)
    test_singele_lset_memory(False)
    test_multiple_lset_memory()
