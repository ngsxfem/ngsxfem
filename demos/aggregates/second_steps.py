from ngsolve import *
from xfem import *
from ngsolve.meshes import MakeStructured2DMesh
SetNumThreads(1)

SetTestoutFile("test_second_step.out")

mesh = MakeStructured2DMesh(quads=False, nx=2, ny=2)

levelset = x - 0.6
gfu = GridFunction(H1(mesh))
InterpolateToP1(levelset, gfu)

ci = CutInfo(mesh, gfu)
roots = ci.GetElementsOfType(NEG)
bads = ci.GetElementsOfType(IF)

EA = ElementAggregation(mesh)
EA.Update(roots, bads)

ba_facets = EA.patch_interior_facets
dw = dFacetPatch(definedonelements=ba_facets)

fes = H1(mesh, order=1)
u, v = fes.TnT()
h = specialcf.mesh_size
bf_dw = 0.1*1.0/h*1.0/h*(u-u.Other())*(v-v.Other()) * dw
lf = v*dx

# with TaskManager():
PatchDummy(EA, fes, fes, bf_dw, lf)

