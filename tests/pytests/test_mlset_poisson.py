# ------------------------------ LOAD LIBRARIES -------------------------------
import pytest
from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.solvers import PreconditionedRichardson as PreRic
from xfem import *
from xfem.mlset import *

ngsglobals.msg_level = 1
SetNumThreads(4)


# -------------------------------- PARAMETERS ---------------------------------
full_order = 4
orders = [1, 2, 3]
h0 = 0.4
mesh_levels = range(0, 6)

rate_tol = 0.125


# ----------------------------------- DATA ------------------------------------
u_ex = 16 * x * (1 - x) * y * (1 - y)
grad_u_ex = CoefficientFunction((16 * (1 - 2 * x) * y * (1 - y),
                                 16 * x * (1 - x) * (1 - 2 * y)))
rhs = 32 * (y * (1 - y) + x * (1 - x))


def level_sets():
    return [-y, x - 1, y - 1, -x]


# --------------------------------- GEOMETRY ----------------------------------
geo = SplineGeometry()
geo.AddRectangle((-0.2, -0.2), (1.2, 1.2),
                 bcs=("bottom", "right", "top", "left"))


# -------------------------------- MAIN RUTINE --------------------------------


def SolvePossionOnUnitSquare(level_sets, mesh, k, rhs, gamma_n=10, gamma_s=0.1, 
        u_ex=None, grad_u_ex=None, inverse="sparsecholesky"):

    nr_ls = len(level_sets)

    # ------------------------- Finite Element Space --------------------------
    V = H1(mesh, order=k, dgjumps=True)

    gfu = GridFunction(V)
    freedofs = BitArray(V.ndof)

    # -------------------------- Levelset & Cut-Info --------------------------
    level_sets_p1 = tuple(GridFunction(H1(mesh, order=1))
                          for i in range(nr_ls))
    for i, lsetp1 in enumerate(level_sets_p1):
        InterpolateToP1(level_sets[i], lsetp1)

    square = DomainTypeArray([(NEG, NEG, NEG, NEG)])
    square.Compress(level_sets_p1)
    lset_dom_inner = {"levelset": level_sets_p1, "domain_type": square}

    boundary = square.Boundary()
    lsets_bnd = {}
    for dtt in boundary:
        lsets_bnd[dtt] = {"levelset": level_sets_p1, "domain_type": dtt}

    mlci = MultiLevelsetCutInfo(mesh, level_sets_p1)

    # ---------------------------- Element Markers ----------------------------
    els_hasneg, els_if = BitArray(mesh.ne), BitArray(mesh.ne)
    els_if_singe = {dtt: BitArray(mesh.ne) for dtt in boundary}
    facets_gp = BitArray(mesh.nedge)

    els_hasneg[:] = False
    els_hasneg |= mlci.GetElementsWithContribution(square)

    els_if[:] = False
    els_if |= els_hasneg & ~mlci.GetElementsOfType(square)

    for i, (dtt, els_bnd) in enumerate(els_if_singe.items()):
        els_bnd[:] = False
        els_bnd |= mlci.GetElementsWithContribution(dtt)

    facets_gp[:] = False
    facets_gp |= GetFacetsWithNeighborTypes(mesh, a=els_hasneg, b=els_if,
                                            use_and=True)

    freedofs[:] = False
    freedofs |= GetDofsOfElements(V, els_hasneg) & V.FreeDofs()

    # --------------------------- (Bi)Linear Forms ----------------------------
    u, v = V.TnT()
    h = specialcf.mesh_size
    normals = square.GetOuterNormals(level_sets_p1)

    diffusion = InnerProduct(Grad(u), Grad(v))

    def nitsche(n):
        return - InnerProduct(Grad(u) * n, v) - InnerProduct(Grad(v) * n, u) \
                    + (gamma_n * k * k / h) * InnerProduct(u, v)
    ghost_penalty = gamma_s / (h**2) * (u - u.Other()) * (v - v.Other())

    forcing = rhs * v

    # ------------------------------ Integrators ------------------------------
    a = RestrictedBilinearForm(V, element_restriction=els_hasneg,
                               facet_restriction=facets_gp, check_unused=False)
    a += SymbolicBFI(lset_dom_inner, form=diffusion,
                     definedonelements=els_hasneg)
    for bnd, n in normals.items():
        a += SymbolicBFI(lsets_bnd[bnd], form=nitsche(n),
                         definedonelements=els_if_singe[bnd])
    a += SymbolicFacetPatchBFI(form=ghost_penalty, skeleton=False,
                               definedonelements=facets_gp)

    f = LinearForm(V)
    f += SymbolicLFI(lset_dom_inner, form=forcing)

    # ----------------------------- Solve Problem -----------------------------
    f.Assemble()
    a.Assemble()
    inv = a.mat.Inverse(freedofs=freedofs,inverse=inverse)

    gfu.vec.data = PreRic(a=a, rhs=f.vec, pre=inv, freedofs=freedofs)

    # --------------------------- Error Computation ---------------------------
    if u_ex:
        err_l2 = sqrt(Integrate(lset_dom_inner, cf=InnerProduct(gfu - u_ex,
                                                                gfu - u_ex),
                                mesh=mesh, order=2 * k))
    else:
        err_l2 = float("NaN")
    if grad_u_ex:
        err_h1 = sqrt(Integrate(lset_dom_inner,
                                cf=InnerProduct(Grad(gfu) - grad_u_ex,
                                                Grad(gfu) - grad_u_ex),
                                mesh=mesh, order=2 * (k - 1)))
    else:
        err_h1 = float("NaN")

    del inv, a, f, gfu
    del els_hasneg, els_if, els_if_singe, facets_gp, freedofs
    del lset_dom_inner, lsets_bnd

    return err_l2, err_h1

# ----------------------------- UTILITY FUNTIONS ------------------------------
def CompLogRate(y1, y2, delta=log(2)):
    return (log(y2) - log(y1)) / delta

# -------------------------------- FIRST CHECK --------------------------------
def test_solution_in_space(): 
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=h0))
        errl2, errh1 = SolvePossionOnUnitSquare(
            level_sets(), mesh, full_order, rhs, u_ex=u_ex, grad_u_ex=grad_u_ex)

        assert errl2 < 1e-12, errh1 < 1e-12

        del mesh


# ----------------------------- CONVERGENCE STUDY -----------------------------
def test_check_convergence_order():
    errs = {"l2": {}, "h1": {}}

    with TaskManager():
        for Lx in mesh_levels:

            ngmesh = geo.GenerateMesh(maxh=h0)
            for i in range(Lx):
                ngmesh.Refine()
            mesh = Mesh(ngmesh)

            for k in orders:
                errl2, errh1 = SolvePossionOnUnitSquare(
                    level_sets(), mesh, k, rhs, u_ex=u_ex, grad_u_ex=grad_u_ex,
                    gamma_s=0.5)

                errs["l2"][(Lx, k)] = errl2
                errs["h1"][(Lx, k)] = errh1

            del ngmesh, mesh


    # ---------------------------- Post-Processing ----------------------------
    for k in orders:
        print("\n Results with k = {:}:".format(k))
        print("--------------------------------------------")
        print("h_max   | L2 err.  | Rate | H1 err.  | Rate ")
        print("--------------------------------------------")

        ratesl2, ratesh1 = 0, 0

        for i, Lx in enumerate(mesh_levels):

            if i > 0:
                rate_l2 = CompLogRate(errs["l2"][(Lx, k)],
                                      errs["l2"][(Lx - 1, k)])
                rate_h1 = CompLogRate(errs["h1"][(Lx, k)],
                                      errs["h1"][(Lx - 1, k)])
                if i > 1:
                    ratesl2 += rate_l2 
                    ratesh1 += rate_h1

                print("{:3.1e}   {:5.3e}  {:4.2f}   {:5.3e}  {:4.2f} ".format(
                    h0 * 0.5**Lx, errs["l2"][(Lx, k)], rate_l2,
                    errs["h1"][(Lx, k)], rate_h1))

            else:
                print("{:3.1e}   {:5.3e}  ----   {:5.3e}  ---- ".format(
                    h0 * 0.5**Lx, errs["l2"][(Lx, k)], errs["h1"][(Lx, k)]))
        print("--------------------------------------------")

        ratesl2 = ratesl2 / (len(mesh_levels) - 2)
        ratesh1 = ratesh1 / (len(mesh_levels) - 2)

        print("Average L2 rate (on final meshes): {}".format(ratesl2))
        print("Average H1 rate (on final meshes): {}".format(ratesh1))

        assert k + 1 - ratesl2 < rate_tol
        assert k - ratesh1 < rate_tol


# -----------------------------------------------------------------------------
# --------------------------- RUN TESTS SEPARATELY ----------------------------
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    test_solution_in_space()
    test_check_convergence_order()
