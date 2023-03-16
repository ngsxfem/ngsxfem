from ngsolve import *
from xfem import *
from ngsolve.meshes import MakeStructured2DMesh
SetNumThreads(1)

SetTestoutFile("test_second_step.out")

mesh = MakeStructured2DMesh(quads=False, nx=3, ny=2)

levelset = x - 0.8
gfu = GridFunction(H1(mesh))
InterpolateToP1(levelset, gfu)

ci = CutInfo(mesh, gfu)
roots = ci.GetElementsOfType(NEG)
bads = ci.GetElementsOfType(IF)

EA = ElementAggregation(mesh)
EA.Update(roots, bads)

ba_facets = EA.patch_interior_facets
dw = dFacetPatch(definedonelements=ba_facets)

#fes = L2(mesh, order=0)
fes = H1(mesh, order=2)
u, v = fes.TnT()
h = specialcf.mesh_size
bf_dw = 0.1*1.0/h*1.0/h*(u-u.Other())*(v-v.Other()) * dw

# with TaskManager():
E = SetupAggEmbedding(EA, fes, bf_dw)

gfu = GridFunction(fes)
Draw(gfu)
baser = E.CreateRowVector()
for i in range(E.width):
    baser[:] = 0
    baser[i] = 1
    gfu.vec.data = E * baser
    Redraw()
    input(i)
#print(E)
