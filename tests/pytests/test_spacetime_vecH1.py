from ngsolve import *
from netgen.geom2d import SplineGeometry
from xfem import *
from math import pi
from xfem.lset_spacetime import *
ngsglobals.msg_level = 1

def test_spacetime_vecH1_diffops():
    SetNumThreads(1)

    # Polynomial order in time
    k_t = 1
    # Polynomial order in space
    k_s = k_t
    # Polynomial order in time for level set approximation
    lset_order_time = k_t
    # Integration order in time
    time_order = 2 * k_t

    maxh = 0.5

    # Outer domain:
    rect = SplineGeometry()
    rect.AddRectangle([-0.5,-0.5], [0.5,0.5])

    # ----------------------------------- MAIN ------------------------------------
    ngmesh = rect.GenerateMesh(maxh=maxh, quad_dominated=False)
    mesh = Mesh(ngmesh)

    # Spatial FESpace for solution
    fes1 = VectorH1(mesh, order=k_s, dgjumps=True, dirichlet=".*")
    # Time finite element (nodal!)
    tfe = ScalarTimeFE(k_t)

    gfs = GridFunction(fes1)

    dst = dxtref(mesh, time_order=time_order)

    st_fes = tfe * fes1

    gfu = GridFunction(st_fes)

    gfs.Set(CF((1 - 3.5 + 3*x + 4*y, 2 - 4 + 5*x + 6*y)))
    gfu.vec.data[0:len(gfs.vec)] = gfs.vec
    gfs.Set(CF((1 + 3.5 + 3*x + 4*y, 2 + 4 + 5*x + 6*y)))
    gfu.vec.data[len(gfs.vec):len(gfu.vec)] = gfs.vec

    print(gfu.vec)

    q = []
    for i, eval in enumerate([
                lambda gfu : gfu[0], 
                lambda gfu : gfu[1],
                lambda gfu : grad(gfu)[0,0],
                lambda gfu : grad(gfu)[0,1],
                lambda gfu : grad(gfu)[1,0],
                lambda gfu : grad(gfu)[1,1],
                lambda gfu : dtref(gfu)[0],
                lambda gfu : dtref(gfu)[1] 
                ]):
        q.append(Integrate( eval(gfu) * dst, mesh))
        assert( abs(i+1-q[i]) < 1e-10 )


    gfs.Set(CF((5*x,4*y)))
    gfu.vec.data[0:len(gfs.vec)] = gfs.vec
    gfs.Set(CF((2+5*x,2+4*y)))
    gfu.vec.data[len(gfs.vec):len(gfu.vec)] = gfs.vec

    q.append(Integrate( div(gfu) * dst, mesh))
    assert( abs(9-q[-1]) < 1e-10 )


    gfs.Set(CF((10 + 3*x + 4*y, 11 + 5*x + 6*y)))
    gfu.vec.data[0:len(gfs.vec)] = gfs.vec
    gfs.Set(CF((12 + 3*x + 4*y, 13 + 5*x + 6*y)))
    gfu.vec.data[len(gfs.vec):len(gfu.vec)] = gfs.vec


    for i, eval in enumerate([
                lambda gfu : fix_tref(gfu,0)[0],
                lambda gfu : fix_tref(gfu,0)[1],
                lambda gfu : fix_tref(gfu,1)[0],
                lambda gfu : fix_tref(gfu,1)[1],
                ]):
        q.append(Integrate( eval(gfu) * dst, mesh))
        assert( abs(i+10-q[-1]) < 1e-10 )


    print(q)