from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh
from numpy import nan 

SetTestoutFile("test.out")

levelset = x - 0.77653
mesh = MakeStructured2DMesh(nx=20, ny=20, quads=True)

EA = ElementAggregation(mesh)

gfu = GridFunction(H1(mesh))
InterpolateToP1(levelset, gfu)

ci = CutInfo(mesh, gfu)
roots = ci.GetElementsOfType(NEG)
bads = ci.GetElementsOfType(IF)
EA.Update(roots, bads)

# print("EA.GetInnerPatchFacets()", EA.GetInnerPatchFacets())
els_surround_patch_facets = GetElementsWithNeighborFacets(mesh, EA.patch_interior_facets)


Draw(BitArrayCF(els_surround_patch_facets), mesh, "surrounding_facets")
Draw(BitArrayCF(EA.els_in_trivial_patch), mesh, "els_in_trivial_patch")
Draw(BitArrayCF(EA.els_in_nontrivial_patch), mesh, "els_in_nontrivial_patch")
Draw(BitArrayCF(roots), mesh, 'roots')
Draw(BitArrayCF(bads), mesh, 'bads')

gf_patch_index = GridFunction(L2(mesh))
gf_patch_index.vec.FV().NumPy()[:] = EA.element_to_patch[:]
#gf_patch_index.vec.FV().NumPy()[gf_patch_index.vec.FV().NumPy()==-1]=nan
Draw(gf_patch_index,mesh,"patchindex")
