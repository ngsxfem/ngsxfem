from ngsolve import *
from xfem import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from sys import argv

import time

ngsglobals.msg_level = 1

def test_fes_timing(dimension=2,stdfes=True,quad=False, order=1):
    if dimension == 2:
        mesh = Mesh (unit_square.GenerateMesh(maxh=0.2,quad_dominated=quad))
    else:
        mesh = Mesh (unit_cube.GenerateMesh(maxh=0.2,quad_dominated=quad))

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
    container.append(("Update",1e9*(te-ts)/steps))
    return container

## dd/mm/yyyy format
date = time.strftime("%Y/%m/%d %H:%M:%S")

if len(argv) > 1:
    basedir = argv[1]
else:
    basedir = "./"
if len(argv) > 2:
    id = str(argv[2])
else:
    id = "0000"

filename1 = basedir + "fes_timings"
for dim in [2,3]:
    filename2 = filename1 + "_dim" + str(dim)
    for stdfes in [True,False]:
        if stdfes:
            filename3 = filename2 + "_std"
        else:
            filename3 = filename2 + "_x"
        for order in [1,3]:
            filename4 = filename3 + "_order" + str(order)
            report = test_fes_timing(dimension=dim,stdfes=stdfes,order=order)
            for key, entry in report:
                filename = filename4 + "_" + key.replace(" ", "_")
                f = open(filename, 'a')
                f.write("{:8}\t{:20}\t{:.6e}\n".format(id,date,entry))
                f.close()
