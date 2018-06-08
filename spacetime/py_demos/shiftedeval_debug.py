# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
from netgen.geom2d import unit_square
from xfem import *

ngsglobals.msg_level = 1

# works for me
mesh = Mesh(unit_square.GenerateMesh(maxh=0.28))

# this crashes
#mesh = Mesh(unit_square.GenerateMesh(maxh=0.29))

# H1-conforming finite element space
fes = H1(mesh, order=3, dirichlet=[1,2,3,4])
fes_dfm = H1(mesh, order=3, dim=2)

gfu_new = GridFunction(fes)

gfu_old = GridFunction(fes)

dfm_back = GridFunction(fes_dfm)
#dfm_back.vec[37:] = 0.1
dfm_back.Set(CoefficientFunction((0.2*sin(5*y),0.2*cos(5*x))))
for i in range(2*mesh.nv):
    dfm_back.vec[i] = 0.0
dfm_forth = GridFunction(fes_dfm)

mesh.SetDeformation(dfm_back)
exact = sin(10*y)
gfu_old.Set(exact)
l2error_old = sqrt (Integrate ( (gfu_old-exact)*(gfu_old-exact), mesh, order=10))
mesh.UnsetDeformation()
Draw (gfu_old,mesh,"gfu_old")
Draw (dfm_back,mesh,"dfmback")
gfu_new.Set(shifted_eval(gfu_old,dfm_back,dfm_forth))
Draw (gfu_new,mesh,"gfu_new")

print ("L2-error(old):", l2error_old)
print ("L2-error(new):", sqrt (Integrate ( (gfu_new-exact)*(gfu_new-exact), mesh, order=10)))
