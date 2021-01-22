from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import sys
import time

ngsglobals.msg_level = 1

def test_fes_timing(dimension=2,stdfes=True,quad_dominated=False, order=1):
    if dimension == 2:
        mesh = Mesh (unit_square.GenerateMesh(maxh=0.2,quad_dominated=quad_dominated))
    else:
        mesh = Mesh (unit_cube.GenerateMesh(maxh=0.2,quad_dominated=quad_dominated))

    Vhs = H1(mesh, order=order, dirichlet=[1,2,3,4])

    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1((sqrt(sqrt(x*x+y*y)) - 1.0),lsetp1)

    if stdfes:
        Vh = Vhs
    else:
        Vhx = XFESpace(Vhs, lsetp1)
        Vh = Vhx
    container = Vh.__timing__()
    ts = time.time()
    steps = 5
    for i in range(steps):
        Vh.Update()
    te = time.time()
    return 1e9*(te-ts)/steps

testcases = [ (2, False, 1, False),
              (2,  True, 1, False),
              (3, False, 1, False),
              (3,  True, 1, False),
              (2, False, 3, False),
              (2,  True, 3, False),
              (3, False, 1, True),
              (3,  True, 1, True),
]

metricname_base = "fes_timings"
for dim, stdfes, order, taskmanager in testcases:
    metricname = metricname_base + "_dim" + str(dim)
    if stdfes:
        metricname += "_std"
    else:
        metricname += "_x"
    metricname += "_order" + str(order)
    if taskmanager:
        metricname += "_TM"
    else:
        metricname += "_NoTM"
    if taskmanager:
        with TaskManager():
            report = test_fes_timing(dimension=dim,stdfes=stdfes,order=order)
    else:
        report = test_fes_timing(dimension=dim,stdfes=stdfes,order=order)
    print(metricname,report, file=sys.stderr)
