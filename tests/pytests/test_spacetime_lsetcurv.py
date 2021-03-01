import pytest
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from math import pi
from xfem.lset_spacetime import *

@pytest.mark.parametrize("imax", [5])
@pytest.mark.parametrize("order", [2,3,4])
def test_spacetime_lsetcurving_maxdist(imax, order):

    ngsglobals.msg_level = 1

    maxdist_at_reflevels = []
    # polynomial order in time
    k_t = order
    # polynomial order in space
    k_s = k_t

    for i in range(imax):
        square = SplineGeometry()
        square.AddRectangle([-0.6,-0.6],[0.6,0.6])
        ngmesh = square.GenerateMesh(maxh=0.5**(i+1), quad_dominated=False)
        mesh = Mesh (ngmesh)

        tstart = 0
        tend = 0.5
        delta_t = (tend - tstart)/(2**(i))

        told = Parameter(tstart)
        t = told + delta_t*tref

        lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = k_t, threshold=0.5, discontinuous_qn=True)

        r0 = 0.5                                 

        rho =  CoefficientFunction((1/(pi))*sin(2*pi*t))
        r = sqrt(x**2+(y-rho)**2)
        levelset= CoefficientFunction(r -r0)

        maxdists = []
        while tend - told.Get() > delta_t/2:
            dfm = lset_adap_st.CalcDeformation(levelset)
            mesh.deformation = dfm
            maxdist = lset_adap_st.CalcMaxDistance(levelset)    
            mesh.deformation = None
            maxdists.append(maxdist)
            told.Set(told.Get() + delta_t)

        print("i: ", i, "\tdist max: ", max(maxdists))
        maxdist_at_reflevels.append(max(maxdists))

    eocs = [log(maxdist_at_reflevels[i-1]/maxdist_at_reflevels[i])/log(2) for i in range(1,len(maxdist_at_reflevels))]
    print("eocs: ", eocs)
    avg_eoc = sum(eocs[-2:])/len(eocs[-2:])
    print("avg: ", avg_eoc)

    assert avg_eoc > order + 0.75

if __name__ == "__main__":
    test_spacetime_lsetcurving_maxdist(5, 2)