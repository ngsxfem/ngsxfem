import xfem
from xfem import check_if_ngsolve_newer_than

def test_ngsolve_version():
    assert(check_if_ngsolve_newer_than(xfem.__ngsolve_required__))