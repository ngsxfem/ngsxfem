from ngsolve import *
from ngsolve.meshes import *
from xfem import *

def test_xfes_ndof_2D():
    mesh = MakeStructured2DMesh(quads=False,nx=4,ny=4,mapping=lambda x,y: (2*x-1,2*y-1))
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1((sqrt(x*x+y*y) - 1.0/3.0),lsetp1)
    Vh = H1(mesh, order=1, dirichlet=[1,2,3,4])
    Vhx = XFESpace(Vh, lsetp1)
    assert Vh.ndof == 25
    assert Vhx.ndof == 7

def test_xfes_ndof_3D():
    mesh = MakeStructured3DMesh(hexes=False,nx=4,ny=4,nz=4,mapping=lambda x,y,z: (2*x-1,2*y-1,2*z-1))
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1((sqrt(x*x+y*y) - 1.0/3.0),lsetp1)
    Vh = H1(mesh, order=1, dirichlet=[1,2,3,4])
    Vhx = XFESpace(Vh, lsetp1)
    assert Vh.ndof == 125
    assert Vhx.ndof == 35
