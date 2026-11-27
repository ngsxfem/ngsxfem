"""
Multigrid for an unfitted interface problem (xfem.cutmg), as in the jupyter tutorial cutmg.ipynb:
P1 CutFEM on a hierarchy of bisected meshes with the prolongation for spaces with inactive dofs.
"""
from math import sqrt as msqrt

from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh
from xfem import *
from xfem.cutmg import LinearMGIterator

mu = [1e-5, 1]
levelset = sqrt(x**2 + y**2) - 1
distsq = x**2 + y**2 - 1
solution = [mu[1] * distsq, mu[0] * distsq]
coef_f = [-4 * mu[0] * mu[1], -4 * mu[0] * mu[1]]


def cut_geometry(mesh):
    lsetp1 = GridFunction(H1(mesh, order=1))
    InterpolateToP1(levelset, lsetp1)
    return CutInfo(mesh, lsetp1), lsetp1


def assemble(mesh, CutVh, ci, lsetp1, gamma_stab=10, lambda_nitsche=10):
    ba_facets = [GetFacetsWithNeighborTypes(mesh, a=ci.GetElementsOfType(dt),
                                            b=ci.GetElementsOfType(IF))
                 for dt in [HASNEG, HASPOS]]
    n_lset = 1.0 / grad(lsetp1).Norm() * grad(lsetp1)
    h = specialcf.mesh_size
    kappa = [mu[1] / (mu[0] + mu[1]), mu[0] / (mu[0] + mu[1])]
    nitsche_penalty = 2 * mu[0] * mu[1] / (mu[0] + mu[1]) * lambda_nitsche
    u, v = CutVh.TnT()
    flux_u = sum(-kappa[i] * mu[i] * grad(u[i]) * n_lset for i in [0, 1])
    flux_v = sum(-kappa[i] * mu[i] * grad(v[i]) * n_lset for i in [0, 1])
    jmp_u, jmp_v = u[0] - u[1], v[0] - v[1]
    doms = [{"levelset": lsetp1, "domain_type": dt} for dt in [NEG, POS, IF]]
    a = BilinearForm(CutVh, symmetric=False)
    for i in [0, 1]:
        a += SymbolicBFI(doms[i], form=mu[i] * grad(u[i]) * grad(v[i]))
        a += SymbolicFacetPatchBFI(form=gamma_stab * mu[i] / h / h
                                   * (u[i] - u[i].Other()) * (v[i] - v[i].Other()),
                                   skeleton=False, definedonelements=ba_facets[i])
    a += SymbolicBFI(doms[2], form=flux_u * jmp_v + flux_v * jmp_u
                     + nitsche_penalty / h * jmp_u * jmp_v)
    a.Assemble()
    f = LinearForm(CutVh)
    for i in [0, 1]:
        f += SymbolicLFI(doms[i], form=coef_f[i] * v[i])
    f.Assemble()
    return a, f


def test_cutmg_interface_poisson():
    nref = 2
    mesh = MakeStructured2DMesh(quads=False, nx=10, ny=10,
                                mapping=lambda x, y: (3 * x - 1.5, 3 * y - 1.5))
    l2errors, mgiter = [], None
    for lvl in range(nref + 1):
        if lvl == 0:
            ci, lsetp1 = cut_geometry(mesh)
            # the spaces are updated by hand after the refinement (active dofs change)
            Vh = H1(mesh, order=1, dirichlet=".*", autoupdate=False)
            VhNeg = Compress(Vh, active_dofs=GetDofsOfElements(Vh, ci.GetElementsOfType(HASNEG)))
            VhPos = Compress(Vh, active_dofs=GetDofsOfElements(Vh, ci.GetElementsOfType(HASPOS)))
            CutVh = FESpace([VhNeg, VhPos], dgjumps=True, autoupdate=False)
            gfu = GridFunction(CutVh)
        else:
            mesh.Refine()
            ci, lsetp1 = cut_geometry(mesh)
            Vh.Update()
            CutVh.components[0].SetActiveDofs(GetDofsOfElements(Vh, ci.GetElementsOfType(HASNEG)))
            CutVh.components[1].SetActiveDofs(GetDofsOfElements(Vh, ci.GetElementsOfType(HASPOS)))
            CutVh.Update()
            gfu.Update()
        a, f = assemble(mesh, CutVh, ci, lsetp1)
        if lvl == 0:
            mgiter = LinearMGIterator(a=a, ci=ci, lsetp1=lsetp1, nref=nref, mesh=mesh, nu=1,
                                      maxit=30, tol=1e-8, printinfo=False,
                                      ProlType=P1Prolongation, ifsolver="cg")
        else:
            mgiter.Update(a, CutVh, ci)
        gfu.components[1].Set(solution[1], BND)
        f.vec.data -= a.mat * gfu.vec
        update = gfu.vec.CreateVector()
        update.data = mgiter * f.vec

        # the multigrid iteration has converged
        proj = Projector(mask=CutVh.FreeDofs(), range=True)
        res = f.vec.CreateVector()
        res.data = proj * (f.vec - a.mat * update)
        rhs = f.vec.CreateVector()
        rhs.data = proj * f.vec
        assert Norm(res) < 1e-6 * Norm(rhs)

        gfu.vec.data += update
        doms = [{"levelset": lsetp1, "domain_type": dt} for dt in [NEG, POS]]
        l2errors.append(msqrt(sum(Integrate(doms[i], cf=(gfu.components[i] - solution[i])**2,
                                            mesh=mesh, order=2) for i in [0, 1])))
    print("L2 errors:", l2errors)
    # second order convergence under refinement
    for coarse, fine in zip(l2errors, l2errors[1:]):
        assert fine < coarse / 3
