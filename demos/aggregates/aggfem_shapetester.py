from ngsolve import *
from xfem import *
from ngsolve.meshes import MakeStructured2DMesh

mesh = MakeStructured2DMesh(quads=False, nx=3, ny=2)

levelset = x - 0.8
gflset = GridFunction(H1(mesh))
InterpolateToP1(levelset, gflset)

ci = CutInfo(mesh, gflset)
roots = ci.GetElementsOfType(NEG)
bads = ci.GetElementsOfType(IF)


#fes = L2(mesh, order=0)
fes = H1(mesh, order=2)

u, v = fes.TnT()
h = specialcf.mesh_size

EA = ElementAggregation(mesh, roots, bads)
E = AggEmbedding(EA, fes)

gfu = GridFunction(fes)
Draw(gfu)
baser = E.CreateRowVector()
for i in range(E.width):
    baser[:] = 0
    baser[i] = 1
    gfu.vec.data = E * baser
    Redraw()
    input(i)
print(" -- done -- ")
