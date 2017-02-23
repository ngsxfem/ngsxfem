from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube

def test_xfes_ndof_2D():
    mesh = Mesh("pytests/mesh2D.vol.gz")
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1((sqrt(x*x+y*y) - 1.0/3.0),lsetp1)
    Vh = H1(mesh, order=1, dirichlet=[1,2,3,4])
    Vhx = XFESpace(Vh, lsetp1)
    assert Vh.ndof == 25
    assert Vhx.ndof == 8

def test_xfes_ndof_3D():
    mesh = Mesh("pytests/mesh3D.vol.gz")
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1((sqrt(x*x+y*y) - 1.0/3.0),lsetp1)
    Vh = H1(mesh, order=1, dirichlet=[1,2,3,4])
    Vhx = XFESpace(Vh, lsetp1)
    assert Vh.ndof == 125
    assert Vhx.ndof == 33
