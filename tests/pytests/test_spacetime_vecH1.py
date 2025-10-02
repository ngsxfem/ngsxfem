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

    #print(gfu.vec)

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

    gfs.Set(CF((14, 15)))
    gfu.vec.data[0:len(gfs.vec)] = gfs.vec
    gfs.Set(CF((18, 19)))
    gfu.vec.data[len(gfs.vec):len(gfu.vec)] = gfs.vec

    for i, t in enumerate([0,0.5,1]):
        for j in range(2):
            q.append(Integrate( fix_tref(gfu,t)[j] * dst, mesh))
            assert( abs(i*2+j+14-q[-1]) < 1e-10 )

    print(q)

import pytest
@pytest.mark.parametrize("k_t", [0,1,2])
@pytest.mark.parametrize("k_s", [1,2,3])
@pytest.mark.parametrize("n_threads", [1,2])
def test_spacetime_vecH1_diffops_compare(k_t, k_s, n_threads):
    solver = "pardiso"
    SetNumThreads(n_threads)
    time_order = 2*k_t+2
    
    from netgen.occ import MoveTo, Circle, Face, OCCGeometry
    base = MoveTo(0,0).Circle(1.1).Face()
    hole = MoveTo(0,0).Circle(1.5/5).Face()
    base -= hole
    geo = OCCGeometry(base, dim=2)
    ngmesh = geo.GenerateMesh(maxh=0.4)
    mesh = Mesh(ngmesh)
    
    mesh.Curve(k_s)

    dst = dxtref(mesh, time_order=time_order)
    V = VectorH1(mesh,order=k_s, dgjumps=True, dirichlet="default")
    tfe = ScalarTimeFE(k_t)

    st_vfes = tfe*V

    u,v = st_vfes.TnT()

    a = BilinearForm(st_vfes)
    a += (grad(u) | grad(v)) * dst 
    a.Assemble()

    W = H1(mesh,order=k_s, dgjumps=True, dirichlet="default")
    st_fes = (tfe*W)*(tfe*W)

    (u1,u2),(v1,v2) = st_fes.TnT()
    us,vs = (u1,u2),(v1,v2)
    
    def sgrad(u):
        return CF((grad(u[0])[0], grad(u[0])[1], grad(u[1])[0], grad(u[1])[1]),dims=(2,2))

    a_s = BilinearForm(st_fes)
    a_s += (sgrad(us) | sgrad(vs)) * dst
    
    a_s.Assemble()

    print(Norm(a.mat.AsVector()))
    print(Norm(a_s.mat.AsVector()))
    assert(abs(Norm(a.mat.AsVector()) - Norm(a.mat.AsVector())) < 1e-10)

    