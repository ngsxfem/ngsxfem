# ngsolve stuff
from ngsolve import *
# For TraceFEM-Integrators (convenience)
from xfem.tracefem import *

def CheckDGpattern(mesh,order,levelset,resultdict):
    """ 
    for a given mesh, a given polynomial degree and a level set function
    determine the sparsity of a discontinuous Galerkin TraceFEM disc. 
    Results are stored in the result dictionary 'resultdict'.
    """
    Vh_l2 = L2(mesh, order=order, dirichlet=[])
    Vh_l2_tr = TraceFESpace(mesh, Vh_l2, levelset, dgjumps=True)
    resultdict["dg_ndofs"] = Vh_l2_tr.ndof
    resultdict["dg_global_ndofs"] = Vh_l2_tr.ndof
    b = BilinearForm(Vh_l2_tr)
    b.Assemble(heapsize=10000000)
    resultdict["dg_nze"] = b.mat.AsVector().size
    resultdict["dg_nze_cond"] = b.mat.AsVector().size

def CheckCGpattern(mesh,order,levelset,resultdict):
    """ 
    for a given mesh, a given polynomial degree and a level set function
    determine the sparsity of a continuous Galerkin TraceFEM disc. 
    (without ghost penalty stabilization terms).
    Results are stored in the result dictionary 'resultdict'.
    """
    Vh_h1 = H1(mesh, order=order, dirichlet=[])
    Vh_h1_tr = TraceFESpace(mesh, Vh_h1, levelset)

    global_cg_ndofs = Vh_h1_tr.ndof
    for i in range(Vh_h1_tr.ndof):
        if (Vh_h1_tr.CouplingType(i) == COUPLING_TYPE.LOCAL_DOF):
            global_cg_ndofs = global_cg_ndofs - 1

    resultdict["cg_ndofs"] = Vh_h1_tr.ndof
    resultdict["cg_global_ndofs"] = global_cg_ndofs

    b = BilinearForm(Vh_h1_tr, flags = {"eliminate_internal" : False})
    b.Assemble(heapsize=10000000)
    resultdict["cg_nze"] = b.mat.AsVector().size
    b = BilinearForm(Vh_h1_tr, flags = {"eliminate_internal" : True})
    b.Assemble(heapsize=10000000)
    resultdict["cg_nze_cond"] = b.mat.AsVector().size

def CheckCGGPpattern(mesh,order,levelset,resultdict):
    """ 
    for a given mesh, a given polynomial degree and a level set function
    determine the sparsity of a continuous Galerkin TraceFEM disc. 
    with ghost penalty stabilization terms.
    Results are stored in the result dictionary 'resultdict'.
    """
    Vh_h1 = H1(mesh, order=order, dirichlet=[])
    Vh_h1_tr = TraceFESpace(mesh, Vh_h1, levelset, dgjumps=True)
    resultdict["cggp_ndofs"] = Vh_h1_tr.ndof
    resultdict["cggp_global_ndofs"] = Vh_h1_tr.ndof
    b = BilinearForm(Vh_h1_tr)
    b.Assemble(heapsize=10000000)
    resultdict["cggp_nze"] = b.mat.AsVector().size
    resultdict["cggp_nze_cond"] = b.mat.AsVector().size

def CheckHDGpattern(mesh,order,levelset,resultdict):
    """ 
    for a given mesh, a given polynomial degree and a level set function
    determine the sparsity of a HDG TraceFEM disc. 
    Results are stored in the result dictionary 'resultdict'.
    """
    Vh_l2 = L2(mesh, order=order, dirichlet=[])
    Vh_l2_tr = TraceFESpace(mesh, Vh_l2, levelset, postpone_update = True)
    Vh_facet = FacetFESpace(mesh, order=order, dirichlet=[], flags = {"highest_order_dc" : False})
    Vh_facet_tr = TraceFESpace(mesh, Vh_facet, levelset, postpone_update = True)
    Vh_tr = FESpace([Vh_l2_tr,Vh_facet_tr])

    global_hdg_ndofs = Vh_tr.ndof
    for i in range(Vh_tr.ndof):
        if (Vh_tr.CouplingType(i) == COUPLING_TYPE.LOCAL_DOF):
            global_hdg_ndofs = global_hdg_ndofs - 1

    resultdict["hdg_ndofs"] = Vh_tr.ndof
    resultdict["hdg_global_ndofs"] = global_hdg_ndofs

    b = BilinearForm(Vh_tr, flags = {"eliminate_internal" : False})
    b.Assemble(heapsize=10000000)
    resultdict["hdg_nze"] = b.mat.AsVector().size

    b = BilinearForm(Vh_tr, flags = {"eliminate_internal" : True})
    b.Assemble(heapsize=10000000)
    resultdict["hdg_nze_cond"] = b.mat.AsVector().size

    Vh_l2 = L2(mesh, order=order, dirichlet=[])
    Vh_l2_tr = TraceFESpace(mesh, Vh_l2, levelset, postpone_update = True)
    Vh_facet = FacetFESpace(mesh, order=order, dirichlet=[], flags = {"highest_order_dc" : True})
    Vh_facet_tr = TraceFESpace(mesh, Vh_facet, levelset, postpone_update = True)
    Vh_tr = FESpace([Vh_l2_tr,Vh_facet_tr])

    global_hdg_ndofs = Vh_tr.ndof
    for i in range(Vh_tr.ndof):
        if (Vh_tr.CouplingType(i) == COUPLING_TYPE.LOCAL_DOF):
            global_hdg_ndofs = global_hdg_ndofs - 1

    resultdict["hdg_hodc_ndofs"] = Vh_tr.ndof
    resultdict["hdg_hodc_global_ndofs"] = global_hdg_ndofs

    b = BilinearForm(Vh_tr, flags = {"eliminate_internal" : False})
    b.Assemble(heapsize=10000000)
    resultdict["hdg_hodc_nze"] = b.mat.AsVector().size

    b = BilinearForm(Vh_tr, flags = {"eliminate_internal" : True})
    b.Assemble(heapsize=10000000)
    resultdict["hdg_hodc_nze_cond"] = b.mat.AsVector().size
