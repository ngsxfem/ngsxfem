# ------------------------------ LOAD LIBRARIES -------------------------------
import pytest
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.lsetcurv import *

ngsglobals.msg_level = 2
SetNumThreads(4)

# -------------------------------- PARAMETERS ---------------------------------
lowerleft, upperright = (0.0,0.0), (2.2,0.41)   # Background domain corners
h_max = 0.07                                    # Mesh diameter
k = 5                                           # Order of velocity space
                                                
mapping = True                                  # Use isoparametric mapping 
nu = 1e-3                                       # Fluid viscosity
gamma_n = 40                                    # Nitsche parameter
gamma_gp = 0.01                                 # Ghost penalty parameter
pReg = 1e-8                                     # Pressure regularisation

inverse = ""                                    # Direct linear solver used
condense = True                                 # Use static condensation 
maxit_newt = 15                                 # Max. nr. of Newton iterations
tol_newt = 1e-10                                # l2 Newton residual tolerance
update_jacobi_tol = 0.15                        # Factor for Jacobean update
reuse_jacobi = True                             # Try to reuse Jacobean
print_newt = True                               # Newton convergence info


# -------------------------- CUTFEM HELPER FUNCTIONS --------------------------
def MarkCutElementsForCondensing(mesh, FES, ghost_facets):
    
    ghost_el = GetElementsWithNeighborFacets(mesh, ghost_facets)

    for el in mesh.Elements():
        if ghost_el[el.nr]:
            for dof_nr in FES.GetDofNrs(el):
                if FES.CouplingType(dof_nr) == COUPLING_TYPE.LOCAL_DOF:
                    FES.SetCouplingType(dof_nr, COUPLING_TYPE.INTERFACE_DOF)

    FES.FinalizeUpdate()


def UpdateMarkers(element_marker, union_elememts, intersection_elements = None):
    
    element_marker.Clear()
    element_marker |= union_elememts
    if intersection_elements:
        element_marker &= intersection_elements

# -------------------------------- NEWTON LOOP --------------------------------
def CutFEM_QuasiNewton(a, alin, u, f, freedofs, maxit=100, maxerr=1e-11,
                       inverse="", dampfactor=1, jacobi_update_tol=0.1,
                       reuse = False, printing=True, **kwargs):

    res = u.vec.CreateVector()
    res_freedofs = u.vec.CreateVector()
    du = u.vec.CreateVector()
    numit = 0
    err,errLast = float("NaN"),float("NaN")
    Updated = "n/a "

    if reuse and not "inv_jacobian" in globals():
        global inv_jacobian
        JacobianAvailable = False
    if reuse and "inv_jacobian" in globals():
        JacobianAvailable = True
    else:
        JacobianAvailable = False

    projector = Projector(freedofs, True)

    if printing: 
        print("Numit\tUpd.J\t||res||_l2")
        print("--------------------------")

    for it in range(maxit):

        a.Assemble()
        res.data = a.mat * u.vec
        if f: res.data -= f
        res_freedofs.data = projector*res
        if it > 0:
            errLast = err
        err = Norm(res_freedofs)
        if printing: print ("{:}\t{:}\t{:1.4e}".format(numit,Updated,err))
        if err < maxerr: break
        elif not JacobianAvailable or err > errLast*jacobi_update_tol:
            UpdateJacobian = True
        else:
            UpdateJacobian = False

        numit += 1

        if UpdateJacobian:
            if JacobianAvailable:
                del inv_jacobian
            alin.Assemble()
            inv_jacobian = alin.mat.Inverse(freedofs, inverse=inverse)
            Updated = True
            JacobianAvailable = True
        else:
            Updated = False

        if alin.condense:
            res.data += alin.harmonic_extension_trans * res
            du.data = inv_jacobian * res
            du.data += alin.inner_solve * res
            du.data += alin.harmonic_extension * du
        else:
            du.data = inv_jacobian*res
            
        u.vec.data -= min(1,numit*dampfactor)*du

    else:
        print("Warning: Newton might not have converged!")
        return (-1,numit)

    if not reuse and JacobianAvailable: del inv_jacobian
    return (0,numit)


def test_st2d1_drag_lift():
    # ----------------------------------- DATA ------------------------------------
    def levelset_func(t):
        return 0.05 - sqrt((x - 0.2) * (x - 0.2) + (y - 0.2) * (y - 0.2))

    u_inflow = CoefficientFunction((4 * 0.3 * y * (0.41 - y) / (0.41**2), 0.0))

    # ------------------------------ BACKGROUND MESH ------------------------------
    geo = SplineGeometry()
    p1, p2, p3, p4, p5, p6 = [geo.AppendPoint(x,y) for x,y in [(0, 0), (0.7, 0), 
                                    (2.2, 0), (2.2, 0.41), (0.7,0.41), (0,0.41)] ]
    geo.Append (["line", p1, p2], leftdomain=1, rightdomain=0, bc="wall")
    geo.Append (["line", p2, p5], leftdomain=1, rightdomain=2)
    geo.Append (["line", p5, p6], leftdomain=1, rightdomain=0, bc="wall")
    geo.Append (["line", p6, p1], leftdomain=1, rightdomain=0, bc="inlet")
    geo.Append (["line", p2, p3], leftdomain=2, rightdomain=0, bc="wall")
    geo.Append (["line", p3, p4], leftdomain=2, rightdomain=0, bc="outlet")
    geo.Append (["line", p4, p5], leftdomain=2, rightdomain=0, bc="wall")

    geo.SetDomainMaxH(1, h_max/6)
    mesh = Mesh(geo.GenerateMesh(maxh=h_max))


    # --------------------------- FINITE ELEMENT SPACE ----------------------------
    V = VectorH1(mesh, order=k, dirichlet="wall|inlet")
    Q = H1(mesh, order=k-1)
    X = FESpace([V,Q], dgjumps=True)

    gfu = GridFunction(X)
    vel,pre = gfu.components
    gfu.vec[:] = 0.0


    # ---------------------------- LEVELSET & CUT-INFO ----------------------------
    # Levelset approximation
    if mapping:
        lset_meshadap = LevelSetMeshAdaptation(mesh, order=k, discontinuous_qn=True)
        deformation = lset_meshadap.CalcDeformation(levelset_func(0.0))
        mesh.SetDeformation(deformation)
        lsetp1 = lset_meshadap.lset_p1
    else:
        lsetp1 = GridFunction(H1(mesh,order=1))                     
        InterpolateToP1(levelset_func(0.0),lsetp1)

    # Integration dictionaries
    lset_neg = {"levelset": lsetp1, "domain_type": NEG, "subdivlvl": 0}
    lset_if = {"levelset": lsetp1, "domain_type": IF, "subdivlvl": 0}

    # Cut-info class
    ci_main = CutInfo(mesh, lsetp1)

    # ------------------------------ ELEMENT MARKERS ------------------------------
    els_hasneg, els_if = BitArray(mesh.ne), BitArray(mesh.ne)
    facets_gp = BitArray(mesh.nedge)
    active_dofs = BitArray(X.ndof)

    UpdateMarkers(els_hasneg, ci_main.GetElementsOfType(HASNEG))
    UpdateMarkers(els_if, ci_main.GetElementsOfType(IF))
    UpdateMarkers(facets_gp, GetFacetsWithNeighborTypes(mesh, a=els_hasneg, 
                                                        b=els_if, use_and=True))

    if condense:
        MarkCutElementsForCondensing(mesh, X, facets_gp)
    UpdateMarkers(active_dofs, GetDofsOfElements(X, els_hasneg), 
                  X.FreeDofs(coupling=condense))

    # ----------------------------- (BI)LINEAR FORMS ------------------------------
    (u, p), (v, q) = X.TnT()

    h = specialcf.mesh_size
    n_levelset = 1.0/Norm(grad(lsetp1)) * grad(lsetp1)

    stokes = nu * InnerProduct(Grad(u), Grad(v)) - p*div(v) - q*div(u)
    convect = InnerProduct(Grad(u) * vel, v)
    convect_lin = InnerProduct(Grad(u) * vel, v) + InnerProduct(Grad(vel) * u, v)

    nitsche  = -nu * InnerProduct(grad(u) * n_levelset,v)
    nitsche += -nu * InnerProduct(grad(v) * n_levelset,u)
    nitsche += nu * (gamma_n * k * k / h) * InnerProduct(u, v)
    nitsche += p * InnerProduct(v, n_levelset)
    nitsche += q * InnerProduct(u, n_levelset)

    ghost_penalty = gamma_gp * nu * (1 / h**2) * (u - u.Other()) * (v - v.Other())
    ghost_penalty += - gamma_gp * (1 / nu) * (p - p.Other()) * (q - q.Other())


    # -------------------------------- INTEGRATORS --------------------------------

    a = RestrictedBilinearForm(X, element_restriction=els_hasneg, 
                               facet_restriction=facets_gp, 
                               check_unused=False, flags={"symmetric": False})
    a += SymbolicBFI(lset_neg, form=stokes + convect, 
                     definedonelements=els_hasneg)
    a += SymbolicBFI(lset_if, form=nitsche, definedonelements=els_if)
    a += SymbolicFacetPatchBFI(form=ghost_penalty, skeleton=False, 
                               definedonelements=facets_gp)

    a_lin = RestrictedBilinearForm(X, element_restriction=els_hasneg, 
                                   facet_restriction=facets_gp, 
                                   check_unused=False, 
                                   flags={"condense":condense, 
                                          "symmetric": False})
    a_lin += SymbolicBFI(lset_neg, form=stokes - pReg * p * q + convect_lin,
                         definedonelements=els_hasneg)
    a_lin += SymbolicBFI(lset_if, form=nitsche, definedonelements=els_if)
    a_lin += SymbolicFacetPatchBFI(form=ghost_penalty, skeleton=False, 
                                   definedonelements=facets_gp)


    # ------------------------- SOLVE STATIONARY PROBLEM --------------------------
    with TaskManager():
        
        vel.Set(u_inflow, definedon=mesh.Boundaries("inlet"))
        CutFEM_QuasiNewton(a=a, alin=a_lin, u=gfu, f=None, freedofs=active_dofs, 
                           maxit=maxit_newt, maxerr=tol_newt, inverse=inverse,
                           jacobi_update_tol=update_jacobi_tol, reuse=reuse_jacobi,
                           printing=print_newt)


    # --------------------------- FUNTIONAL EVALUATION ----------------------------

    drag_x_test, drag_y_test = GridFunction(X), GridFunction(X)
    drag_x_test.components[0].Set(CoefficientFunction((1.0,0.0)))
    drag_y_test.components[0].Set(CoefficientFunction((0.0,1.0)))

    n = specialcf.normal(mesh.dim)
    a_test = BilinearForm(X, symmetric=False, check_unused=False)
    a_test += SymbolicBFI(lset_neg, form=convect, definedonelements=els_hasneg)
    a_test += SymbolicBFI(form=-InnerProduct(nu*Grad(u)*n - p*n,v), skeleton=True, 
                          definedon=mesh.Boundaries("inlet|wall"))
    with TaskManager():
        a_test.Assemble()

    Um = 2 * 0.3 / 3
    C_drag = -2.0 / (0.1 * Um**2) * a_test(gfu, drag_x_test)
    C_lift = -2.0 / (0.1 * Um**2) * a_test(gfu, drag_y_test)
    pdiff = pre(0.15, 0.2) - pre(0.25, 0.2)

    err_drag = abs(C_drag - 5.57953523384)
    err_lift = abs(C_lift - 0.010618948146)
    err_p = abs(pdiff - 0.11752016697)

    print("\n     C_drag      C_lift     pdiff")
    print("Val: {:10.8f}  {:10.8f}  {:10.8f}".format(C_drag, C_lift, pdiff))
    print("Err: {:1.2e}    {:1.2e}    {:1.2e}".format(err_drag, err_lift, err_p))

    assert err_drag < 2.5e-5
    assert err_lift < 2e-6
    assert err_p    < 2e-4

if __name__ == "__main__":
    test_st2d1_drag_lift()
