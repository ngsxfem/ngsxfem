from netgen.geom2d import unit_square
from ngsolve import *
from xfem import *

mesh = Mesh(unit_square.GenerateMesh(maxh=0.5))
fes = H1(mesh, order=2)

EA = ElementAggregation(mesh)
els_all, els_none = BitArray(mesh.ne), BitArray(mesh.ne)
els_all.Set()
els_none.Clear() 
EA.Update(els_all, els_none)

print('non trivial patch elements\n', EA.els_in_nontrivial_patch)
print('elements in trivial patches\n', EA.els_in_trivial_patch)

gfu = GridFunction(fes)
gfu.Set(x**3 + sin(y))

u,v = fes.TnT()
bf = u*v*dx
lf = gfu*v*dx

gfupatch = GridFunction(fes)
gfupatch.vec.data = PatchwiseSolve(EA,fes,bf,lf)

Draw(gfu,mesh,"gfuset")
Draw(gfupatch,mesh,"gfupatch")
print(sqrt(Integrate(InnerProduct(gfu-gfupatch,gfu-gfupatch) * dx ,mesh)))
