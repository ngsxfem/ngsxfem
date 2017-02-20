from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube

def test_xfes_ndof_2D(quad=False):
    mesh = Mesh (unit_square.GenerateMesh(maxh=2,quad_dominated=quad))
    mesh.Refine()
    mesh.Refine()
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1((sqrt(x*x+y*y) - 1.0/3.0),lsetp1)
    Vh = H1(mesh, order=1, dirichlet=[1,2,3,4])
    Vhx = XFESpace(Vh, mesh, lsetp1)
    assert Vh.ndof == 65
    assert Vhx.ndof == 10

def test_xfes_ndof_3D(quad=False):
    mesh = Mesh (unit_cube.GenerateMesh(maxh=2,quad_dominated=quad))
    mesh.Refine()
    mesh.Refine()
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1((sqrt(x*x+y*y) - 1.0/3.0),lsetp1)
    Vh = H1(mesh, order=1, dirichlet=[1,2,3,4])
    Vhx = XFESpace(Vh, mesh, lsetp1)
    assert Vh.ndof == 125
    assert Vhx.ndof == 33
