# CutFEM with MPI (without facet terms): the active dofs of Compress and Restrict are
# consistent across the ranks, the operator and the solution do not depend on the
# partition. Every rank also computes the sequential reference.
# Run with  mpirun -np N python3 -m pytest test_cutfem_mpi.py  (N = 2, 3, 4, ...)
import pickle

import pytest
mpi = pytest.importorskip("mpi4py.MPI")
import netgen.meshing
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve import *
from pyngcore import MPI_Comm
from xfem import *
from xfem.lsetcurv import LevelSetMeshAdaptation


def maxh(dim):
    return 0.1 if dim == 2 else 0.25


def meshes(dim):
    comm = MPI_Comm(mpi.COMM_WORLD)
    geo = unit_square if dim == 2 else unit_cube
    ngmesh = geo.GenerateMesh(maxh=maxh(dim)) if comm.rank == 0 else None
    seq = Mesh(mpi.COMM_WORLD.bcast(ngmesh))
    if comm.size == 1:
        par = Mesh(ngmesh)
    elif comm.rank == 0:
        par = Mesh(pickle.loads(pickle.dumps(ngmesh)).Distribute(comm))
    else:
        par = Mesh(netgen.meshing.Mesh.Receive(comm))
    return seq, par


def quantities(mesh, order, deform, restrict):
    dim = mesh.dim
    levelset = sqrt((x - 0.5)**2 + (y - 0.5)**2
                    + ((z - 0.5)**2 if dim == 3 else 0)) - 0.31
    if deform:
        lsetadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000,
                                          discontinuous_qn=True)
        deformation = lsetadap.CalcDeformation(levelset)
        lsetp1 = lsetadap.lset_p1
    else:
        deformation = None
        lsetp1 = GridFunction(H1(mesh, order=1))
        InterpolateToP1(levelset, lsetp1)

    ci = CutInfo(mesh, lsetp1)
    hasneg = ci.GetElementsOfType(HASNEG)
    hasif = ci.GetElementsOfType(IF)

    Vhbase = H1(mesh, order=order)
    if restrict:
        Vh = Restrict(Vhbase, hasneg)
    else:
        Vh = Compress(Vhbase, GetDofsOfElements(Vhbase, hasneg))
    u, v = Vh.TnT()
    h = specialcf.mesh_size
    nl = Normalize(grad(lsetp1))
    dxc = dCut(lsetp1, NEG, definedonelements=hasneg, deformation=deformation)
    dsc = dCut(lsetp1, IF, definedonelements=hasif, deformation=deformation)

    lam = 10 * order**2
    a = BilinearForm(Vh, symmetric=False)
    a += grad(u) * grad(v) * dxc
    a += (-grad(u) * nl * v - grad(v) * nl * u + lam / h * u * v) * dsc
    # instead of a ghost penalty (which needs facet terms across the ranks)
    a += 1e-3 * u * v * dx(definedonelements=hasneg)
    f = LinearForm(Vh)
    f += v * dxc
    a.Assemble()
    f.Assemble()

    gfu = GridFunction(Vh)
    if mesh.comm.size > 1:
        inv = a.mat.Inverse(Vh.FreeDofs(), inverse="masterinverse")
    else:
        inv = a.mat.Inverse(Vh.FreeDofs())
    gfu.vec.data = inv * f.vec
    w = GridFunction(Vh)
    w.Set(x * x + y)
    tmp = w.vec.CreateVector()
    tmp.data = a.mat * w.vec
    return [Vh.ndofglobal, InnerProduct(tmp, w.vec), InnerProduct(gfu.vec, f.vec)]


@pytest.mark.parametrize("dim,deform", [(2, False), (2, True), (3, False)])
@pytest.mark.parametrize("restrict", [False, True])
def test_cutfem_partition_independent(dim, deform, restrict):
    order = 2
    seq, par = (quantities(m, order, deform, restrict) for m in meshes(dim))
    assert seq[0] == par[0], (seq, par)
    for s, p in zip(seq[1:], par[1:]):
        assert abs(s - p) <= 1e-10 * abs(s), (seq, par)
