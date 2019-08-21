# integration on lset domains
from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

from ngsolve.meshes import *

import pytest

@pytest.mark.parametrize("quad", [False,True])
@pytest.mark.parametrize("order", [1,2,3])

def test_intcurved(quad, order):
    levelset = sqrt(x*x+y*y)-0.5
    referencevals = { POS : 4.0-0.25*pi, NEG : 0.25*pi, IF : pi }

    N=4
    errors_uncurved = dict()
    errors_curved = dict()
    eoc_uncurved = dict()
    eoc_curved = dict()

    for key in [NEG,POS,IF]:
        errors_curved[key] = []
        errors_uncurved[key] = []
        eoc_curved[key] = []
        eoc_uncurved[key] = []

    refinements = 5
    if order == 1:
        refinements += 2
    for reflevel in range(refinements):

        mesh = MakeStructured2DMesh(quads = quad, nx=N*2**reflevel, ny=N*2**reflevel,
                                    mapping = lambda x,y : (2*x-1,2*y-1))    
        lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)
        lsetp1 = lsetmeshadap.lset_p1

        f = CoefficientFunction (1.0)

        for key in [NEG,POS,IF]:
            # Applying the mesh deformation
            deformation = lsetmeshadap.CalcDeformation(levelset)

            integrals_uncurved = Integrate(levelset_domain = { "levelset" : lsetp1, "domain_type" : key},
                                           cf=f, mesh=mesh, order = order)

            mesh.SetDeformation(deformation)
            integrals_curved = Integrate(levelset_domain = { "levelset" : lsetp1, "domain_type" : key},
                                           cf=f, mesh=mesh, order = order)
            # Unapply the mesh deformation (for refinement)
            mesh.UnsetDeformation()

            errors_curved[key].append(abs(integrals_curved - referencevals[key]))
            errors_uncurved[key].append(abs(integrals_uncurved - referencevals[key]))
        # refine cut elements:
        # if not quad:
        #     RefineAtLevelSet(gf=lsetmeshadap.lset_p1)

    for key in [NEG,POS,IF]:
        eoc_curved[key] = [log(a/b)/log(2) for (a,b) in zip (errors_curved[key][0:-1],errors_curved[key][1:]) ]
        eoc_uncurved[key] = [log(a/b)/log(2) for (a,b) in zip (errors_uncurved[key][0:-1],errors_uncurved[key][1:]) ]

    # print("errors (uncurved):  \n{}\n".format(errors_uncurved))
    # print("   eoc (uncurved):  \n{}\n".format(   eoc_uncurved))
    print("errors (  curved):  \n{}\n".format(  errors_curved))
    print("   eoc (  curved):  \n{}\n".format(     eoc_curved))

    print("avg.eoc(curved, IF):  \n{}\n".format(     sum(eoc_curved[IF][2:])/len(eoc_curved[IF][2:])))
    print("avg.eoc(curved,NEG):  \n{}\n".format(     sum(eoc_curved[NEG][2:])/len(eoc_curved[NEG][2:])))
    print("avg.eoc(curved,POS):  \n{}\n".format(     sum(eoc_curved[POS][2:])/len(eoc_curved[POS][2:])))

    if (order > 1):
        assert errors_curved[IF][-1] < 1e-5
        assert errors_curved[NEG][-1] < 1e-5
        assert errors_curved[POS][-1] < 1e-5
    else:
        assert errors_curved[IF][-1] < 1e-4
        assert errors_curved[NEG][-1] < 1e-4
        assert errors_curved[POS][-1] < 1e-4

    s = 0
    if order == 1:
        s+=2
    assert sum(eoc_curved[IF][s:])/len(eoc_curved[IF][s:]) > order + 0.75
    assert sum(eoc_curved[NEG][s:])/len(eoc_curved[NEG][s:]) > order + 0.75
    assert sum(eoc_curved[POS][s:])/len(eoc_curved[POS][s:]) > order + 0.75
