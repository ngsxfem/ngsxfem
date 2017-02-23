from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube

ngsglobals.msg_level = 1

def test_xfes_timing_2D(quad=False):
    mesh = Mesh (unit_square.GenerateMesh(maxh=0.2,quad_dominated=quad))
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1((sqrt(sqrt(x*x+y*y)) - 1.0),lsetp1)
    Vh = H1(mesh, order=1, dirichlet=[1,2,3,4])
    Vhx = XFESpace(Vh, lsetp1)
    return Vhx.__timing__()

def test_xfes_timing_3D(quad=False):
    mesh = Mesh (unit_cube.GenerateMesh(maxh=0.2,quad_dominated=quad))
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1((sqrt(sqrt(x*x+y*y)) - 1.0),lsetp1)
    Vh = H1(mesh, order=1, dirichlet=[1,2,3,4])
    Vhx = XFESpace(Vh, lsetp1)
    return Vhx.__timing__()

print("\nCalculations: \n\n")
timing2d = test_xfes_timing_2D()
timing3d = test_xfes_timing_3D()

print("\n\nReport:")
print("XFES Timings 2D:\n",timing2d)
print("XFES Timings 3D:\n",timing3d)
