from ngsolve import *
from xfem import *
from netgen.geom2d import SplineGeometry
from netgen.geom2d import unit_square

import pytest

@pytest.mark.parametrize("space", [H1,L2])
@pytest.mark.parametrize("order", [1,2,3])

def test_patch_local_solve(space,order):
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.5))
    fes = space(mesh, order=order)

    EA = ElementAggregation(mesh)
    els_all, els_none = BitArray(mesh.ne), BitArray(mesh.ne)
    els_all.Set()
    els_none.Clear() 
    EA.Update(els_all, els_none)

    #print('non trivial patch elements\n', EA.els_in_nontrivial_patch)
    #print('elements in trivial patches\n', EA.els_in_trivial_patch)

    gfu = GridFunction(fes)
    gfu.Set(x**3 + sin(y))

    u,v = fes.TnT()
    bf = u*v*dx
    lf = gfu*v*dx

    gfupatch = GridFunction(fes)
    gfupatch.vec.data = PatchwiseSolve(EA,fes,bf,lf)

    Draw(gfu,mesh,"gfuset")
    Draw(gfupatch,mesh,"gfupatch")
    diff = sqrt(Integrate(InnerProduct(gfu-gfupatch,gfu-gfupatch) * dx ,mesh))
    assert(diff < 1e-12)


@pytest.mark.parametrize("order", [1, 2])
def test_patch_local_solve_facet_restriction(order):
    """The definedonelements of a facet patch integral refer to facets. On the
    patches, only the ghost penalty determines the dofs of the cut elements:
    the mass term lives only on the root elements (definedonelements=NEG).
    The patchwise solution is then the polynomial that was projected."""
    from ngsolve.meshes import MakeStructured2DMesh
    mesh = MakeStructured2DMesh(quads=False, nx=8, ny=8)
    lsetp1 = GridFunction(H1(mesh, order=1))
    InterpolateToP1(x - 0.51, lsetp1)
    ci = CutInfo(mesh, lsetp1)
    EA = ElementAggregation(mesh)
    EA.Update(ci.GetElementsOfType(NEG), ci.GetElementsOfType(IF))
    assert sum(EA.patch_interior_facets) > 0

    fes = L2(mesh, order=order, dgjumps=True)
    u, v = fes.TnT()
    # polynomial of degree <= order: exactly representable in the space
    poly = 1 + x + 2 * y + (x * y - x**2 if order > 1 else 0)
    dxinner = dx(definedonelements=ci.GetElementsOfType(NEG))
    dP = dFacetPatch(definedonelements=EA.patch_interior_facets)
    bf = u * v * dxinner + (u - u.Other()) * (v - v.Other()) * dP
    lf = poly * v * dxinner

    gfu = GridFunction(fes)
    gfu.vec.data = PatchwiseSolve(EA, fes, bf, lf)
    dxactive = dx(definedonelements=ci.GetElementsOfType(HASNEG))
    err = sqrt(Integrate((gfu - poly)**2 * dxactive, mesh))
    assert err < 1e-9
