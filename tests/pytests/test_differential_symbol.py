import pytest
from netgen.geom2d import unit_square
from ngsolve import *
from xfem import *
from xfem.lsetcurv import *
ngsglobals.msg_level = 0

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
levelset = (x - 0.5)**2 + (y - 0.5)**2 - 0.33**2


@pytest.mark.parametrize('order', [1, 2, 3])
@pytest.mark.parametrize('DOM', [POS, IF, NEG])
@pytest.mark.parametrize('subdivlvl', [0, 1, 2])
def test_cut_symbols_straightcut(order, DOM, subdivlvl):

    if order > 1 and subdivlvl == 0:
        lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order,
                                              levelset=levelset)
        lset_h = lsetmeshadap.lset_p1
    elif order == 1 and subdivlvl == 0:
        lsetmeshadap = NoDeformation(mesh, levelset)
        lset_h = lsetmeshadap.lset_p1
    elif subdivlvl > 0:
        lsetmeshadap = NoDeformation(mesh, levelset)
        lset_h = GridFunction(H1(mesh, order=order))
        lset_h.Set(levelset)

    deform = lsetmeshadap.deform

    ci = CutInfo(mesh, lset_h)
    els_hasdom = ci.GetElementsOfType(HAS(DOM))
    facets_none = BitArray(mesh.nedge)
    facets_none.Clear()

    V = H1(mesh, order=order)
    u, v = V.TnT()

    gfu = GridFunction(V)
    gfu.Set(sin(x))

    w1, w2 = gfu.vec.CreateVector(), gfu.vec.CreateVector()

    # Bilinear Form
    form = Grad(u) * Grad(v) + u * v

    # Differential Symbol Version
    a1 = RestrictedBilinearForm(V, element_restriction=els_hasdom,
                                facet_restriction=facets_none,
                                check_unused=False)
    a1 += form * dCut(lset_h, DOM, subdivlvl=subdivlvl,
                      definedonelements=els_hasdom, deformation=deform)
    a1.Assemble()
    w1.data = a1.mat * gfu.vec

    # SymbolicBFI version
    lset_dom = {'levelset': lset_h, 'domain_type': DOM, 'subdivlvl': subdivlvl}
    a2 = RestrictedBilinearForm(V, element_restriction=els_hasdom,
                                facet_restriction=facets_none,
                                check_unused=False)
    a2 += SymbolicBFI(levelset_domain=lset_dom, form=form,
                      definedonelements=els_hasdom)
    mesh.SetDeformation(deform)
    a2.Assemble()
    mesh.UnsetDeformation()
    w2.data = a2.mat * gfu.vec

    # Check
    w1.data -= w2
    diff = Norm(w1)
    print(f'diff : {diff}')
    assert diff < 1e-12


@pytest.mark.parametrize('order', [1, 2, 3])
@pytest.mark.parametrize('DOM', [POS, NEG])
@pytest.mark.parametrize('skeleton', [True, False])
@pytest.mark.parametrize('element_boundary', [True, False])
def test_cut_symbols_dg(order, DOM, skeleton, element_boundary):
    if skeleton is False and element_boundary is False:
        return None
    if order > 1:
        lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order,
                                              levelset=levelset)
    else:
        lsetmeshadap = NoDeformation(mesh, levelset)

    lset_h = lsetmeshadap.lset_p1
    deform = lsetmeshadap.deform

    ci = CutInfo(mesh, lset_h)
    els_hasdom = ci.GetElementsOfType(HAS(DOM))
    facets_dom = GetFacetsWithNeighborTypes(mesh, a=els_hasdom, b=els_hasdom)
    els_none = BitArray(mesh.ne)
    els_none[:] = False

    V = H1(mesh, order=order, dgjumps=True)
    u, v = V.TnT()

    gfu = GridFunction(V)
    gfu.Set(sin(x))

    w1, w2 = gfu.vec.CreateVector(), gfu.vec.CreateVector()

    # Bilinear form to integrate
    nF = specialcf.normal(mesh.dim)
    flux_u = -0.5 * (grad(u) + grad(u.Other())) * nF
    flux_v = -0.5 * (grad(v) + grad(v.Other())) * nF
    jump_u = u - u.Other()
    jump_v = v - v.Other()

    form = jump_u * jump_v + flux_u * jump_v + flux_v * jump_u

    # Differential Symbol Version
    a1 = RestrictedBilinearForm(V, element_restriction=els_hasdom,
                                facet_restriction=facets_dom,
                                check_unused=False)
    a1 += form * dCut(lset_h, DOM, definedonelements=facets_dom,
                      deformation=deform, skeleton=skeleton,
                      element_boundary=element_boundary)
    a1.Assemble()
    w1.data = a1.mat * gfu.vec

    # SymbolicBFI version
    lset_dom = {'levelset': lset_h, 'domain_type': DOM, 'subdivlvl': 0}
    a2 = RestrictedBilinearForm(V, element_restriction=els_hasdom,
                                facet_restriction=facets_dom,
                                check_unused=False)
    a2 += SymbolicBFI(levelset_domain=lset_dom, form=form,
                      definedonelements=facets_dom, skeleton=skeleton,
                      element_boundary=element_boundary)
    mesh.SetDeformation(deform)
    a2.Assemble()
    mesh.UnsetDeformation()
    w2.data = a2.mat * gfu.vec

    # Check
    w1.data -= w2
    diff = Norm(w1)
    print(f'diff : {diff}')
    assert diff < 1e-12


@pytest.mark.parametrize('order', [1, 2, 3])
@pytest.mark.parametrize('DOM', [POS, NEG])
def test_cut_symbols_facetpatch(order, DOM):

    # Higher order level set approximation
    if order > 1:
        lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order,
                                              levelset=levelset)
    else:
        lsetmeshadap = NoDeformation(mesh, levelset)

    deform = lsetmeshadap.deform
    lsetp1 = lsetmeshadap.lset_p1

    ci = CutInfo(mesh, lsetp1)
    els_hasdom = ci.GetElementsOfType(HAS(DOM))
    els_if = ci.GetElementsOfType(IF)
    ba_facets = GetFacetsWithNeighborTypes(mesh, a=els_hasdom, b=els_if)
    els_gp = GetElementsWithNeighborFacets(mesh, ba_facets)

    V = H1(mesh, order=order, dgjumps=True)
    u, v = V.TnT()

    gfu = GridFunction(V)
    gfu.Set(sin(x))

    w1, w2 = gfu.vec.CreateVector(), gfu.vec.CreateVector()

    # Bilinear form
    form = (u - u.Other()) * (v - v.Other())

    # Differential Symbol Version
    a1 = RestrictedBilinearForm(V, element_restriction=els_gp,
                                facet_restriction=ba_facets,
                                check_unused=False)
    a1 += form * dFacetPatch(definedonelements=ba_facets, deformation=deform)
    a1.Assemble()
    w1.data = a1.mat * gfu.vec

    # SymbolicBFI version
    a2 = RestrictedBilinearForm(V, element_restriction=els_gp,
                                facet_restriction=ba_facets,
                                check_unused=False)
    a2 += SymbolicFacetPatchBFI(form=form, skeleton=False,
                                definedonelements=ba_facets)
    mesh.SetDeformation(deform)
    a2.Assemble()
    mesh.UnsetDeformation()
    w2.data = a2.mat * gfu.vec

    # Check
    w1.data -= w2
    diff = Norm(w1)
    print(f'diff : {diff}')
    assert diff < 1e-12


@pytest.mark.parametrize('order', [1, 2, 3])
@pytest.mark.parametrize('DOM', [POS, IF, NEG])
@pytest.mark.parametrize('subdivlvl', [0, 1, 2])
def test_cut_symbol_integrate(order, DOM, subdivlvl):

    if order > 1 and subdivlvl == 0:
        lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order,
                                              levelset=levelset)
        lset_h = lsetmeshadap.lset_p1
    elif order == 1 and subdivlvl == 0:
        lsetmeshadap = NoDeformation(mesh, levelset)
        lset_h = lsetmeshadap.lset_p1
    elif subdivlvl > 0:
        lsetmeshadap = NoDeformation(mesh, levelset)
        lset_h = GridFunction(H1(mesh, order=order))
        lset_h.Set(levelset)

    deform = lsetmeshadap.deform

    f = x**2 * sin(y)

    # Differential Symbol version
    dx = dCut(levelset=lset_h, domain_type=DOM, order=order,
              subdivlvl=subdivlvl, deformation=deform)
    integral_dx = Integrate(f * dx, mesh=mesh)

    # Old version
    lset_dom = {'levelset': lset_h, 'domain_type': DOM, 'subdivlvl': subdivlvl}
    mesh.SetDeformation(deform)
    integral = Integrate(levelset_domain=lset_dom, cf=f, mesh=mesh,
                         order=order)
    mesh.UnsetDeformation()

    assert abs(integral_dx - integral) < 1e-12
