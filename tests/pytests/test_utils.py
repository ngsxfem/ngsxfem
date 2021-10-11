from test_spacetimecutrule import test_spacetime_spaceP4_timeDGP4
from xfem import ngsxfemglobals

def test_multiply_all_eps():
    e1 = test_spacetime_spaceP4_timeDGP4()

    ngsxfemglobals.MultiplyAllEps(1e4)

    e2 = test_spacetime_spaceP4_timeDGP4()

    assert( abs(e1-e2) > 1e-10)
