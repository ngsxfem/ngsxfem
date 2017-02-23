from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube

import time

ngsglobals.msg_level = 1

def test_fes_timing_2D(quad=False, order=1):
    mesh = Mesh (unit_square.GenerateMesh(maxh=0.2,quad_dominated=quad))
    Vh = H1(mesh, order=order, dirichlet=[1,2,3,4])
    container = Vh.__timing__()
    ts = time.time()
    steps = 5
    for i in range(steps):
        Vh.Update()
    te = time.time()
    container.append(("Update",1e9*(te-ts)/steps))
    return container

def test_xfes_timing_2D(quad=False, order=1):
    mesh = Mesh (unit_square.GenerateMesh(maxh=0.2,quad_dominated=quad))
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1((sqrt(sqrt(x*x+y*y)) - 1.0),lsetp1)
    Vh = H1(mesh, order=order, dirichlet=[1,2,3,4])
    Vhx = XFESpace(Vh, lsetp1)
    container = Vhx.__timing__()
    ts = time.time()
    steps = 5
    for i in range(steps):
        Vhx.Update()
    te = time.time()
    container.append(("Update",1e9*(te-ts)/steps))
    return container

def test_fes_timing_3D(quad=False, order=1):
    mesh = Mesh (unit_cube.GenerateMesh(maxh=0.2,quad_dominated=quad))
    Vh = H1(mesh, order=order, dirichlet=[1,2,3,4])
    container = Vh.__timing__()
    ts = time.time()
    steps = 5
    for i in range(steps):
        Vh.Update()
    te = time.time()
    container.append(("Update",1e9*(te-ts)/steps))
    return container

def test_xfes_timing_3D(quad=False, order=1):
    mesh = Mesh (unit_cube.GenerateMesh(maxh=0.2,quad_dominated=quad))
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1((sqrt(sqrt(x*x+y*y)) - 1.0),lsetp1)
    Vh = H1(mesh, order=order, dirichlet=[1,2,3,4])
    Vhx = XFESpace(Vh, lsetp1)
    container = Vhx.__timing__()
    ts = time.time()
    steps = 5
    for i in range(steps):
        Vhx.Update()
    te = time.time()
    container.append(("Update",1e9*(te-ts)/steps))
    return container

print("\nCalculations: \n\n")
stiming2d_k1 = test_fes_timing_2D()
xtiming2d_k1 = test_xfes_timing_2D()
stiming3d_k1 = test_fes_timing_3D()
xtiming3d_k1 = test_xfes_timing_3D()
stiming2d_k3 = test_fes_timing_2D(order=3)
xtiming2d_k3 = test_xfes_timing_2D(order=3)
stiming3d_k3 = test_fes_timing_3D(order=3)
xtiming3d_k3 = test_xfes_timing_3D(order=3)
with TaskManager():
    stiming3d_k1_tm = test_fes_timing_3D()
    xtiming3d_k1_tm = test_xfes_timing_3D()

print("\n\nReport:")
print(" FES Timings 2D - k = 1 :\n",stiming2d_k1)
print("XFES Timings 2D - k = 1 :\n",xtiming2d_k1)
print(" FES Timings 3D - k = 1 :\n",stiming3d_k1)
print("XFES Timings 3D - k = 1 :\n",xtiming3d_k1)
print(" FES Timings 3D - k = 1 with TaskManager :\n",stiming3d_k1_tm)
print("XFES Timings 3D - k = 1 with TaskManager :\n",xtiming3d_k1_tm)
print(" FES Timings 2D - k = 3 :\n",stiming2d_k3)
print("XFES Timings 2D - k = 3 :\n",xtiming2d_k3)
print(" FES Timings 3D - k = 3 :\n",stiming3d_k3)
print("XFES Timings 3D - k = 3 :\n",xtiming3d_k3)
