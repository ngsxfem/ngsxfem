import pytest
from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh


@pytest.mark.parametrize("struc_mesh", [True, False])
@pytest.mark.parametrize("quad", [False])
@pytest.mark.parametrize("ROOTS", [POS, NEG])
@pytest.mark.parametrize("dim, levelset", [(2, x - 0.77654),
                                           (2, (x - 0.5)**2 + (y - 0.5)**2 - 0.3**2),
                                           (3, (x - 0.5)**2 + (y - 0.5)**2 + (z - 0.5)**2 - 0.3**2)
                                           ])
def test_aggregates(dim, struc_mesh, quad, ROOTS, levelset):
    if dim == 2:
        if struc_mesh:
            mesh = MakeStructured2DMesh(nx=20, ny=20, quads=quad)
        else:
            mesh = Mesh(unit_square.GenerateMesh(maxh=0.05,
                                                 quad_dominated=quad))
    if dim == 3:
        if struc_mesh:
            mesh = MakeStructured3DMesh(hexes=quad, nx=10)
        else:
            if quad:
                return None
            cube = CSGeometry()
            cube.Add(OrthoBrick(Pnt(-1, -1, -1), Pnt(1, 1, 1)))
            mesh = Mesh(cube.GenerateMesh(maxh=0.1, quad_dominated=quad))
    EA = ElementAggregation(mesh)

    gfu = GridFunction(H1(mesh))
    InterpolateToP1(levelset, gfu)

    ci = CutInfo(mesh, gfu)
    roots = ci.GetElementsOfType(ROOTS)
    bads = ci.GetElementsOfType(IF)

    EA.Update(roots, bads)

    # print("EA.GetInnerPatchFacets()", EA.GetInnerPatchFacets())
    els_surround_patch = GetElementsWithNeighborFacets(mesh, EA.patch_interior_facets)

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

    return None


@pytest.mark.parametrize("struc_mesh", [True, False])
@pytest.mark.parametrize("quad", [True, False])
@pytest.mark.parametrize("nr", [10])
@pytest.mark.parametrize("dim, levelset", [(2, x - 0.77654),
                                           (2, (x - 0.5)**2 + (y - 0.5)**2 - 0.3**2),
                                           (3, (x - 0.5)**2 + (y - 0.5)**2 + (z - 0.5)**2 - 0.3**2)
                                           ])
def test_elements_with_contribution(dim, struc_mesh, quad, nr, levelset):
    if dim == 2:
        if struc_mesh:
            mesh = MakeStructured2DMesh(nx=20, ny=20, quads=quad)
        else:
            mesh = Mesh(unit_square.GenerateMesh(maxh=0.05,
                                                 quad_dominated=quad))
    if dim == 3:
        if struc_mesh:
            mesh = MakeStructured3DMesh(hexes=quad, nx=10)
        else:
            if quad:
                return None
            cube = CSGeometry()
            cube.Add(OrthoBrick(Pnt(-1, -1, -1), Pnt(1, 1, 1)))
            mesh = Mesh(cube.GenerateMesh(maxh=0.1, quad_dominated=quad))

    gfu = GridFunction(H1(mesh))
    InterpolateToP1(levelset, gfu)
    ci = CutInfo(mesh, gfu)

    for i in range(1, nr):
        els_test = ci.GetElementsWithThresholdContribution(NEG, i / nr)
        assert sum(els_test & ci.GetElementsOfType(POS)) == 0

        els_test = ci.GetElementsWithThresholdContribution(POS, i / nr)
        assert sum(els_test & ci.GetElementsOfType(NEG)) == 0

    return None

if __name__ == "__main__":
    test_aggregates(2, False, False, NEG, (x - 0.5)**2 + (y-0.5)**2 - 0.3**2)
    test_aggregates(3, False, False, NEG, (3, x - 0.77653))

