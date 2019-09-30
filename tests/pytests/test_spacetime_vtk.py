import pytest
from math import pi
from ngsolve import *
from netgen.geom2d import SplineGeometry
from xfem import *

@pytest.mark.parametrize("quad", [True, False])
def test_spacetime_vtk(quad):
  tref = ReferenceTimeVariable()
  square = SplineGeometry()
  A = 1.25
  square.AddRectangle([-A,-A],[A,A])
  ngmesh = square.GenerateMesh(maxh=0.3, quad_dominated=False)
  mesh = Mesh (ngmesh)
  
  #### expression for the time variable: 
  
  coef_told = Parameter(0)
  coef_delta_t = Parameter(0)
  tref = ReferenceTimeVariable()
  t = coef_told + coef_delta_t*tref
  
  #### the data: 
  # radius of disk (the geometry)
  r0 = 0.5
  x0 = lambda t: r0 * cos(pi*t)
  y0 = lambda t: r0 * sin(pi*t)
  
  #levelset= x - t
  levelset= sqrt((x-x0(t))**2+(y-y0(t))**2) - 0.4
  # spatial FESpace for solution
  fesp1 = H1(mesh, order=1)
  # polynomial order in time for level set approximation
  lset_order_time = 1
  # integration order in time
  time_order = 2
  # time finite element (nodal!)
  tfe = ScalarTimeFE(1) 
  # space-time finite element space
  st_fes = SpaceTimeFESpace(fesp1,tfe)
  
  coef_delta_t.Set(0.25)
  
  lset_p1 = GridFunction(st_fes)
  SpaceTimeInterpolateToP1(levelset,tref,lset_p1)
  
  vtk = SpaceTimeVTKOutput(ma=mesh,coefs=[levelset,lset_p1],
                           names=["levelset","lset_p1"],
                           filename="spacetime_vtk_test",
                           subdivision_x=1,subdivision_t=1)
  vtk.Do(t_start=coef_told.Get(), t_end=coef_told.Get() + coef_delta_t.Get())
  import os
  os.remove("spacetime_vtk_test_0.vtk")

if __name__ == "__main__":
    test_spacetime_vtk(False)
    test_spacetime_vtk(True)
    
