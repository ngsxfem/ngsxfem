import pytest
from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square

#@pytest.mark.parametrize("quad", [False])
def test_aggregates():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.1, quad_dominated=False))
    EA = ElementAggregation(mesh)

    gfu = GridFunction(H1(mesh))

    levelset = (x - 0.77654)

    gfu.Set(levelset)
    ci = CutInfo(mesh, gfu)
    roots = ci.GetElementsOfType(NEG)
    bads = ci.GetElementsOfType(IF)

    EA.Update(roots, bads)
    # Draw(gfu)

    print("EA.GetInnerPatchFacets()", EA.GetInnerPatchFacets())
    patch_facets = EA.GetInnerPatchFacets()

    ba_surround_facets = GetElementsWithNeighborFacets(mesh, patch_facets)
    Draw(BitArrayCF(ba_surround_facets), mesh, "surrounding_facets")

    assert sum(ba_surround_facets & ci.GetElementsOfType(POS)) == 0
    facets_gp = GetFacetsWithNeighborTypes(mesh, a=ci.GetElementsOfType(HASNEG),
                                           b=ci.GetElementsOfType(IF), use_and=True, bnd_val_a=False, bnd_val_b=False)
    gp_surround = GetElementsWithNeighborFacets(mesh, facets_gp)

    Draw(BitArrayCF(gp_surround), mesh, 'surround_gp')
    Draw(BitArrayCF(bads), mesh, 'if')

    Draw(BitArrayCF(ba_surround_facets & ~gp_surround), mesh, 'bad_els')

    assert sum(ba_surround_facets & ~gp_surround) == 0


if __name__ == "__main__":
    SetTestoutFile("out.out")
    test_aggregates()
