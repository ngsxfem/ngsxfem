import pytest
from math import pi
from ngsolve import *
from netgen.geom2d import unit_square
from xfem import *

@pytest.mark.parametrize("order", [0,1,2,3,4,5])
def test_spacetime_set(order):
    ngmesh = unit_square.GenerateMesh(maxh=0.8, quad_dominated=False)
    mesh = Mesh (ngmesh)

    coef_told = Parameter(0)
    delta_t = 1
    tref = ReferenceTimeVariable()
    t = coef_told + delta_t*tref

    k_s = k_t = order
    fes1 = H1(mesh, order=k_s)

    tfe = ScalarTimeFE(k_t) 
    tfe_i = ScalarTimeFE(k_t, skip_first_node=True) # interior
    tfe_e = ScalarTimeFE(k_t, only_first_node=True) # exterior (inital values)

    st_fes = SpaceTimeFESpace(fes1,tfe, flags = {"dgjumps": True})
    st_fes_i = SpaceTimeFESpace(fes1,tfe_i, flags = {"dgjumps": True})
    st_fes_e = SpaceTimeFESpace(fes1,tfe_e, flags = {"dgjumps": True})

    gfu_slice = GridFunction(fes1)

    gfu = GridFunction(st_fes)
    gfu_i = GridFunction(st_fes_i)
    gfu_e = GridFunction(st_fes_e)

    for gf in [gfu,gfu_i,gfu_e]:
        V = gf.space
        gf.Set( 1 + t  )
        for i,that in enumerate(V.TimeFE_nodes()):
            if V.IsTimeNodeActive(i):
                val = 1 + that
            else:
                val = 0
            RestrictGFInTime(gf,that,gfu_slice)
            avg = Integrate(gfu_slice,mesh)
            assert(abs(avg-val) < 1e-12)

if __name__ == "__main__":
    test_spacetime_set(2)
    
