from ngsolve import TimeFunction

@TimeFunction
def FindElementClusters(mesh,inner_el,cut_el,report=False):
    # collect initial elements that are associated to clusters (will grow)
    facets = GetFacetsWithNeighborTypes(mesh,a=inner_el,b=cut_el)
    cluster_el = BitArray(mesh.ne)
    cluster_el = inner_el & GetElementsWithNeighborFacets(mesh,facets)

    Draw(BitArrayCF(inner_el),mesh,"inner_el")
    Draw(BitArrayCF(cut_el),mesh,"cut_el")
    Draw(BitArrayCF(cluster_el),mesh,"cluster_el")

    # element-to-cluster association
    el_to_cluster = [None for el in range(mesh.ne)]
    nclusters = 0
    for elnr in range(mesh.ne):
        if cluster_el[elnr]:
            el_to_cluster[elnr] = nclusters
            nclusters += 1

    # cluster-to-master association
    cluster_to_els = [None for c in range(nclusters)]
    cluster = 0
    for elnr in range(mesh.ne):
        if cluster_el[elnr]:
            cluster_to_els[cluster] = [elnr]
            cluster += 1
    

    # elements that are not yet associated to clusters (but will)
    need_cluster_el = BitArray(mesh.ne)
    need_cluster_el.Clear()
    need_cluster_el |= cut_el


    its = 0
    while (sum(need_cluster_el)) > 0:
        if report:
            print("iteration", its)
            print("number of elements that need a cluster:", sum(need_cluster_el))
        next_cluster = []
        for el in mesh.Elements():
            if cluster_el[el.nr]:
                for f in el.facets:
                    for nel in mesh[f].elements:
                        if nel.nr != el.nr and need_cluster_el[nel.nr]:
                            next_cluster.append((nel.nr,el_to_cluster[el.nr]))
        # print(next_cluster)
        for elnr, cluster in next_cluster:
            need_cluster_el[elnr] = False
            if (el_to_cluster[elnr] == None):  ### <- design choice
                el_to_cluster[elnr] = cluster
            cluster_el[elnr] = True
        its += 1

    for elnr in range(mesh.ne):
        cluster = el_to_cluster[elnr]
        if cluster != None:
            if cluster_to_els[cluster][0] != elnr: # <- not the master element 
                cluster_to_els[cluster].append(elnr)

    ##report:
    #number of clusters:
    if report:
        print("number of clusters:", nclusters)
        for cluster_index, els in enumerate(cluster_to_els):
            print("cluster", cluster_index, "contains elements ", els)


        gfci = GridFunction(L2(mesh))
        gfci.vec[:] = float("nan")
        for cluster_index, els in enumerate(cluster_to_els):
            gfc = GridFunction(L2(mesh))
            gfc.vec[:] = float("nan")
            for el in els:
                gfci.vec[el] = cluster_index
                gfc.vec[el] = 1
            Draw(gfc,mesh,"cluster"+str(cluster_index))
        Draw(gfci,mesh,"cluster_index")
            
    return cluster_to_els

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la
from numpy.linalg import inv 

# R = sp.lil_matrix((n_activedofs,Vh.ndof))
# for i,j in enumerate(activedof_to_dof):
#     R[i,j]=1
# #print(R)    
# Ar = R*A*R.transpose()
# Fr = R*F

def SetupClusterMatrix(els_in_cluster, Vh):
    cluster_el = BitArray(mesh.ne); cluster_el.Clear()
    for elnr in els_in_cluster:
        cluster_el[elnr] = True
        
    ba_facets = GetFacetsWithNeighborTypes(mesh,a=cluster_el,b=cluster_el)
    Vhext = FESpace(Vh.type, Vh.mesh, order=Vh.globalorder, dgjumps=True)    
    u,v = Vhext.TnT()
    h = specialcf.mesh_size
    aext = BilinearForm(Vhext, symmetric=False, check_unused=False) # 
    aext += 0.1*1.0/h*1.0/h*(u-u.Other())*(v-v.Other()) * dFacetPatch(definedonelements=ba_facets)
    # aext += SymbolicFacetPatchBFI(form=0.1*1.0/h*1.0/h*(u-u.Other())*(v-v.Other()), skeleton=False, definedonelements=ba_facets)
    aext.Assemble()
    rows,cols,vals = aext.mat.COO()
    S = sp.csr_matrix((vals,(rows,cols)))
    # if sum(cluster_el) > 1:
    #     print(sum(cluster_el))
    #     print("cluster_el: ", cluster_el)
    #     dofs = GetDofsOfElements(Vhext, cluster_el)
    #     print("Array dofs: ", dofs)
    #     print(S[dofs,:][:, dofs].toarray())
    return S

# def CalcDofOverlapping(clusters,Vh):
#     gfmult = GridFunction(Vh)
#     gfmult.vec[:] = 0
    
#     for cluster in clusters:
#         dof_in_cluster = BitArray(Vh.ndof); dof_in_cluster.Clear()
#         for elnr in cluster:
#             for d in Vh.GetDofNrs(ElementId(VOL,elnr)):
#                 dof_in_cluster[d] = True
#         for i, isin in enumerate(dof_in_cluster):
#             if isin:
#                 gfmult.vec[i] += 1
#     print(gfmult.vec)
#     return gfmult.vec


def SetupProlongationMatrix(inner_el, clusters, Vh, report = False):
    
    mesh = Vh.mesh

    active_dof = BitArray(Vh.ndof)
    active_dof.Clear()
    
    for el in Vh.Elements():
        if inner_el[el.nr]:
            for d in el.dofs:
                active_dof[d] = True

    n_active_dof = sum(active_dof)
    print(f'number of active dofs = {n_active_dof}')
    dof2rdof = [-1 for i in range(Vh.ndof)]
    rdof2dof = [-1 for i in range(n_active_dof)]

    n_active_dof = 0
    for i in range(Vh.ndof):
        if active_dof[i]:
            dof2rdof[i] = n_active_dof
            rdof2dof[n_active_dof] = i
            n_active_dof += 1
    if report:
        print("dof2rdof:",dof2rdof)
        print("rdof2dof:",rdof2dof)
                
    P = sp.lil_matrix((Vh.ndof,n_active_dof))

    
    for rdof,dof in enumerate(rdof2dof):
        P[dof,rdof] = 1
    
    for cluster_index, els_in_cluster in enumerate(clusters):
        inner_dofs = []
        ext_dofs = []
        print("cluster", cluster_index)
        for el_in_cluster, elnr in enumerate(els_in_cluster):
            if el_in_cluster == 0: #master element
                for d in Vh.GetDofNrs(ElementId(VOL,elnr)):
                    if d not in inner_dofs:
                        inner_dofs.append(d)
            else: #dependent element
                for d in Vh.GetDofNrs(ElementId(VOL,elnr)):
                    if d not in inner_dofs:
                        if d not in ext_dofs:
                            ext_dofs.append(d)
        # print(inner_dofs)
        # print(ext_dofs)
        S = SetupClusterMatrix(els_in_cluster, Vh)
        Soo = S[ext_dofs,:][:,ext_dofs].todense()
        Soi = S[ext_dofs,:][:,inner_dofs].todense()
        # want Soi * ui + Soo * uo = 0 
        Q =  - inv(Soo) * Soi
        # uo = Q * ui
        # print(Q)

        for i,edof in enumerate(ext_dofs):
            for j,idof in enumerate(inner_dofs):
                P[edof,dof2rdof[idof]] = Q[i,j]
        #break
        
    # print(P)
    # print(P.todense())
    return P



def TestProlongationMatrix(P, inner_el, clusters, Vh, report = False):
    mesh = Vh.mesh
    active_dof = BitArray(Vh.ndof)
    active_dof.Clear()
    
    for el in Vh.Elements():
        if inner_el[el.nr]:
            for d in el.dofs:
                active_dof[d] = True

    n_active_dof = sum(active_dof)
    dof2rdof = [-1 for i in range(Vh.ndof)]
    rdof2dof = [-1 for i in range(n_active_dof)]

    n_active_dof = 0
    for i in range(Vh.ndof):
        if active_dof[i]:
            dof2rdof[i] = n_active_dof
            rdof2dof[n_active_dof] = i
            n_active_dof += 1
    if report:
        print("dof2rdof:",dof2rdof)
        print("rdof2dof:",rdof2dof)
    
    for cluster_index, els_in_cluster in enumerate(clusters):
        dofs = []
        print("cluster", cluster_index)
        for el_in_cluster, elnr in enumerate(els_in_cluster):
            for d in Vh.GetDofNrs(ElementId(VOL,elnr)):
                if d not in dofs:
                    dofs.append(d)
                            
        # print(dofs)


        S = SetupClusterMatrix(els_in_cluster, Vh)
        SS = S[dofs,:][:,dofs]
        for i,dof in enumerate(dofs):
            if dof2rdof[dof] > 0:
                # Q = P.getcol(dof2rdof[dof]).todense().transpose())
                Q = P[dofs,dof2rdof[dof]].toarray().reshape((len(dofs),1))
                n = np.linalg.norm(SS*Q)
                if n > 1e-10:
                    if report:
                        print("test failed with ... ")
                        print("cluster_index",cluster_index)
                        print("dof2rdof[dof]",dof2rdof[dof])
                        print("dofs",dofs)
                    return False
    return True

from ngsolve import *
from xfem import *
if __name__ == "__main__":
    # from netgen import gui
    from ngsolve.meshes import MakeStructured2DMesh
    mesh = MakeStructured2DMesh(quads=False,nx=2,ny=2)
    
    # levelset = y*y+x*x-0.5
    levelset = x-0.6
    # levelset = x+y-1.1
    
    lsetp1 = GridFunction(H1(mesh))
    lsetp1.Set(levelset)
    ci = CutInfo(mesh,lsetp1)
    clusters = FindElementClusters(mesh, ci.GetElementsOfType(NEG), ci.GetElementsOfType(IF), report=False)
    print(clusters)

    #Vh = HDiv(mesh, order=0)
    # Vh = L2(mesh, order=1)
    Vh = H1(mesh, order=1)
    P = SetupProlongationMatrix(ci.GetElementsOfType(NEG), clusters, Vh, report = True)
    print(P)
    # test = TestProlongationMatrix(P, ci.GetElementsOfType(NEG), clusters, Vh, report = False)
    # if not test:
    #     raise Exception("test failed")

    ### Visualize prolongated/agglomerated basis:
    # active_dof = GetDofsOfElements(Vh,ci.GetElementsOfType(NEG))
    # gf_ext_basis = GridFunction(Vh)
    # Draw(gf_ext_basis,mesh,"gf_ext_basis")
    # Draw(gf_ext_basis.Deriv(),mesh,"gf_ext_basis.Deriv()")
    # active_dof_n = 0
    # for dof in range(Vh.ndof):
    #     if active_dof[dof]:
    #         gf_ext_basis.vec[:] = 0
    #         gf_ext_basis.vec.FV().NumPy()[:] = P.getcol(active_dof_n).todense().transpose()
    #         Redraw()
    #         input("active dof number: " + str(active_dof_n))            
    #         active_dof_n += 1
