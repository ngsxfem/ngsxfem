import xfem
from xfem.ngs_check import check_if_ngsolve_newer_than, __ngsolve_required__

def test_ngsolve_version():
    assert(check_if_ngsolve_newer_than(xfem.__ngsolve_required__))