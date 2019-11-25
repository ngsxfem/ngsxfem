# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

import pytest
from ngsolve import *
from ngsolve.meshes import *
from netgen.geom2d import unit_square
from xfem import *

ngsglobals.msg_level = 1

@pytest.mark.parametrize("dimension", [2,3])
def test_shifteval(dimension):
  if dimension == 2:
    mesh = MakeStructured2DMesh(quads = False, nx=8, ny=8)
  else:
    mesh = MakeStructured3DMesh(hexes = False, nx=8, ny=8, nz=8)
    
  # H1-conforming finite element space
  fes = H1(mesh, order=3, dirichlet=[1,2,3,4])
  fes_dfm = H1(mesh, order=3, dim=dimension)
  
  gfu_new = GridFunction(fes)
  
  gfu_old = GridFunction(fes)
  
  dfm_back = GridFunction(fes_dfm)
  #dfm_back.vec[37:] = 0.1
  if dimension == 2:
    dfm_back.Set(CoefficientFunction((0.2*sin(5*y),0.2*cos(5*x))))
  else:
    dfm_back.Set(CoefficientFunction((0.15*sin(5*y),0.15*cos(5*z),0.15*sin(5*x))))
    
  for i in range(dimension*mesh.nv):
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

  error_old = l2error_old
  error_new = sqrt (Integrate ( (gfu_new-exact)*(gfu_new-exact), mesh, order=10))
  print ("L2-error(old):", error_old)
  print ("L2-error(new):", error_new)
  assert error_old < 1e-3
  assert error_new < 1e-3
