import pytest
from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh


# @pytest.mark.parametrize("dim", [2, 3])
# @pytest.mark.parametrize("struc_mesh", [True, False])
# @pytest.mark.parametrize("quad", [True, False])
# @pytest.mark.parametrize("ROOTS", [POS, NEG])
# @pytest.mark.parametrize("levelset" [x - 0.77654, (x - 0.5)**2 + (y-0.5)**4 - 0.3**3])
def test_aggregates(dim, struc_mesh, quad, ROOTS, levelset):
    if dim == 2:
        if struc_mesh:
            mesh = MakeStructured2DMesh(nx=20, ny=20, quads=quad)
        else:
            mesh = Mesh(unit_square.GenerateMesh(maxh=0.05,
                                                 quad_dominated=quad))
    if dim == 3:
        if struc_mesh:
            mesh = MakeStructured2DMesh(nx=20, )
        else:
            mesh = 
    EA = ElementAggregation(mesh)

    gfu = GridFunction(H1(mesh))
    InterpolateToP1(levelset, gfu)

    ci = CutInfo(mesh, gfu)
    roots = ci.GetElementsOfType(ROOTS)
    bads = ci.GetElementsOfType(IF)

    EA.Update(roots, bads)

    # print("EA.GetInnerPatchFacets()", EA.GetInnerPatchFacets())
    patch_facets = EA.GetInnerPatchFacets()
    els_surround_patch = GetElementsWithNeighborFacets(mesh, patch_facets)

    facets_gp = GetFacetsWithNeighborTypes(mesh,
                                           a=ci.GetElementsOfType(HAS(ROOTS)),
                                           b=ci.GetElementsOfType(IF),
                                           use_and=True, bnd_val_a=False,
                                           bnd_val_b=False)
    els_surround_gp = GetElementsWithNeighborFacets(mesh, facets_gp)

    Draw(BitArrayCF(els_surround_patch), mesh, "surrounding_facets")
    Draw(BitArrayCF(els_surround_gp), mesh, 'surround_gp')
    Draw(BitArrayCF(bads), mesh, 'if')
    Draw(BitArrayCF(els_surround_patch & ~els_surround_gp), mesh, 'bad_els')

    assert sum(els_surround_patch & ~els_surround_gp) == 0
    assert sum(els_surround_patch & ~ci.GetElementsOfType(HAS(ROOTS))) == 0

if __name__ == "__main__":
    test_aggregates(2, False, False, NEG, (x - 0.5)**2 + (y-0.5)**2 - 0.3**2)

