# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
from netgen.geom2d import unit_square
from xfem import *

ngsglobals.msg_level = 1

# generate a triangular mesh of mesh-size 0.2
mesh = Mesh(unit_square.GenerateMesh(maxh=0.5))

# H1-conforming finite element space
fes = H1(mesh, order=5, dirichlet=[1,2,3,4])
fes_dfm = H1(mesh, order=3, dim=2)

gfu = GridFunction(fes)

gfu1 = GridFunction(fes)

dfm_back = GridFunction(fes_dfm)
#dfm_back.vec[37:] = 0.1
dfm_back.Set(CoefficientFunction((0.05*sin(5*y),0.05*cos(5*x))))
for i in range(2*mesh.nv):
    dfm_back.vec[i] = 0.0
dfm_forth = GridFunction(fes_dfm)

mesh.SetDeformation(dfm_back)
gfu.Set(sin(10*y))
mesh.UnsetDeformation()
gfu1.Set(shifted_eval(gfu,dfm_back,dfm_forth))
Draw (gfu,mesh,"gfu")
Draw (gfu1,mesh,"gfu1")
Draw (dfm_back,mesh,"dfmback")

#gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
# print (u.vec)

#
## plot the solution (netgen-gui only)
#Draw (gfu)
#Draw (-grad(gfu), mesh, "Flux")
#
#exact = 16*x*(1-x)*y*(1-y)
#print ("L2-error:", sqrt (Integrate ( (gfu-exact)*(gfu-exact), mesh)))
