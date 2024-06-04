/// from ngsxfem
#include "../xfem/cutinfo.hpp"
#include "../xfem/aggregates.hpp"


#include <unordered_map>

using namespace ngsolve;
using namespace ngfem;

namespace ngcomp
{

  ElementAggregation::ElementAggregation (shared_ptr<MeshAccess> ama)
    : ma(ama)
  {
    //std::cout << "Hello from element aggregation" << endl;
  }

  void ElementAggregation::Update(shared_ptr<BitArray> & root_els, shared_ptr<BitArray> & bad_els, LocalHeap & lh)
  {
    const size_t ne = ma->GetNE();
    const size_t nf = ma->GetNFacets();

    el_roots = root_els;

    if (ma->GetCommunicator().Size() > 1)
      throw Exception("ElementAggregation::Update:: No Aggregation implementation for MPI yet");

    /*
    BitArray fine_facet(nf);
    fine_facet.Clear();
    IterateRange
      (ma->GetNE(VOL), lh,
       [&] (int elnr, LocalHeap & lh)
       {
         Array<int> fanums(0,lh);
         fanums = ma->GetElFacets (ElementId(VOL,elnr));
         for (int j=0; j<fanums.Size(); j++)
           fine_facet.SetBitAtomic(fanums[j]);
       });
    */


    shared_ptr<BitArray> facets = GetFacetsWithNeighborTypes(ma,root_els,bad_els,false,false,true,lh);
    shared_ptr<BitArray>patch_root_elements = make_shared<BitArray>(ne);
    patch_root_elements->Clear();
    //cout << *patch_root_elements << endl;
    *patch_root_elements |= *root_els;
    //cout << *patch_root_elements << endl;
    *patch_root_elements &= *GetElementsWithNeighborFacets(ma,facets,lh);
    //cout << *patch_root_elements << endl;
    size_t nc = 0;
    for (size_t i : Range(ne)){
      if (patch_root_elements->Test(i)){
        nc ++;
      }
    }
    (*testout) << "Number of (nontrivial) patches: " << nc << endl;

    element_to_patch.SetSize(ne);
    element_to_patch = -1;

    facet_to_patch.SetSize(nf);
    facet_to_patch = -1;

    patch_roots.SetSize(nc);
    size_t counter = 0;
    for (size_t i: Range(ma->GetNE())){
      if (patch_root_elements->Test(i)){
        element_to_patch[i] = counter;
        patch_roots[counter]=i;
        counter++;
      }
      else
        element_to_patch[i] = -1;
    }

    BitArray front_elements(ne);
    BitArray new_front_elements(ne);
    new_front_elements = *patch_root_elements;

    BitArray not_covered_yet(ne);
    not_covered_yet.Clear();
    not_covered_yet |= *bad_els;
    BitArray newly_covered(ne);
    bool uncovered_elements_left = true;

    inner_patch_facets = make_shared<BitArray> (nf);
    inner_patch_facets->Clear();

    while (not_covered_yet.NumSet() > 0)
    {
      front_elements = new_front_elements;
      new_front_elements.Clear();
      newly_covered.Clear();
      //IterateRange(ne, lh, [&] (int elnr, LocalHeap & lh)
      for (size_t elnr = 0; elnr < ne; elnr++)
      {
        size_t patch_id = element_to_patch[elnr];
        if (front_elements.Test(elnr))
        {
          Array<int> fanums(0,lh);
          fanums = ma->GetElFacets (ElementId(VOL,elnr));
          for (int facnr : fanums)
          {
            Array<int> elnums(0,lh);
            ma->GetFacetElements (facnr, elnums);
            for (int nelnr : elnums)
            {
              if (nelnr == elnr)
                continue; // only consider the neighbor elements
              if (not_covered_yet.Test(nelnr) && (!newly_covered.Test(nelnr)))
              {
                inner_patch_facets->SetBitAtomic(facnr);
                facet_to_patch[facnr]=patch_id;
                new_front_elements.SetBitAtomic(nelnr);
                element_to_patch[nelnr] = patch_id;
                newly_covered.SetBitAtomic(nelnr);
              }
              else
              {
                if (element_to_patch[elnr] == element_to_patch[nelnr])
                {
                  inner_patch_facets->SetBitAtomic(facnr);
                  facet_to_patch[facnr]=patch_id;
                }
              }
            }
          }

        }
      }
      not_covered_yet &= ~newly_covered;
    }

    (*testout) << "element_to_patch = \n" << element_to_patch << endl;

    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto elnr : Range(ne))
      {
        if (element_to_patch[elnr] >= 0 && patch_roots[element_to_patch[elnr]] != elnr)
          creator.Add (element_to_patch[elnr], elnr);
      }
    patch_leafs = creator.MoveTable();

    TableCreator<size_t> fcreator;
    for ( ; !fcreator.Done(); fcreator++)
      for (auto facnr : Range(nf))
        {
          if (facet_to_patch[facnr] >= 0)
            fcreator.Add (facet_to_patch[facnr], facnr);
        }
    patch_facets = fcreator.MoveTable();

    TableCreator<size_t> pcreator;
    for ( ; !pcreator.Done(); pcreator++)
      for (auto pnr : Range(nc))
        {
          pcreator.Add(pnr, patch_roots[pnr]);
          for ( size_t i : patch_leafs[pnr])
            pcreator.Add(pnr, i);
        }
    patch = pcreator.MoveTable();

    (*testout) << "patch id to included element: \n " << patch << endl;
    (*testout) << "patch id to supporting element: \n " << patch_roots << endl;
    (*testout) << "patch id to supported elements: \n " << patch_leafs << endl;
    (*testout) << "patch id to supported connecting facets: \n " << patch_facets << endl;
    (*testout) << "inner patch facets: \n " << *inner_patch_facets << endl;

    n_els_in_nontrivial_patches = 0;
    n_nontrivial_patches = nc;

    n_trivial_els = 0;
    for (size_t i: Range(ma->GetNE())){
      if (root_els->Test(i))
        n_trivial_els ++;
    }
    n_trivial_els -= n_nontrivial_patches;

    el_in_nontrivial_patch = make_shared<BitArray>(ne);
    el_in_nontrivial_patch->Clear();
    for (size_t i: Range(nc))
    {
      n_els_in_nontrivial_patches ++;
      el_in_nontrivial_patch->SetBitAtomic(patch_roots[i]);
      for (size_t j : patch_leafs[i])
      {
        n_els_in_nontrivial_patches ++;
        el_in_nontrivial_patch->SetBitAtomic(j);
      }
    }

    el_in_trivial_patch = make_shared<BitArray>(ne);
    el_in_trivial_patch->Clear();

    *el_in_trivial_patch |= *root_els;
    *el_in_trivial_patch &= ~(*el_in_nontrivial_patch);

    
    trivial_patch_to_element.SetSize(n_trivial_els);
    int cnt = 0;
    for (int i : Range(ne))
      if (el_in_trivial_patch->Test(i))
        trivial_patch_to_element[cnt++] = i;   

  }

  void ElementAggregation::GetPatch(int patchnr, Array<size_t> & ret)
  {
    const int triv_patch_idx = patchnr - GetNNontrivialPatches();
    if (triv_patch_idx < 0)
    {
      ret.SetSize(patch[patchnr].Size());
      ret = patch[patchnr];
    }
    else
    {
      ret.SetSize(1);
      ret = trivial_patch_to_element[triv_patch_idx];
    }
  }

  void ElementAggregation::GetInnerPatchFacets(int patchnr, Array<size_t> & ret)
  {
    const int triv_patch_idx = patchnr - GetNNontrivialPatches();
    if (triv_patch_idx < 0)
    {
      ret.SetSize(patch_facets[patchnr].Size());
      ret = patch_facets[patchnr];
    }
    else
    {
      ret.SetSize0();
    }
  }

  template <typename TFUNC>
  void PatchLoop (shared_ptr<ElementAggregation> elagg, 
                  bool include_trivial_patches,
                  shared_ptr<FESpace> fes_trial,
                  shared_ptr<FESpace> fes_test,
                  shared_ptr<SumOfIntegrals> bf,
                  shared_ptr<SumOfIntegrals> lf,
                  LocalHeap & clh,
                  const TFUNC & patchwise_func
                  )
  {
    (*testout) << "Hello from PatchLoop " << endl;
    size_t npatch = elagg->GetNPatches(include_trivial_patches);
    const BitArray & freedofs_trial = *fes_trial->GetFreeDofs();
    if (freedofs_trial.NumSet() < fes_trial->GetNDof())
      throw Exception("cannot handle non-trivial freedofs array");

    shared_ptr<MeshAccess> ma = elagg->GetMesh();
    Array<size_t> ret;
    for (int i : Range(npatch) )
    {
      elagg->GetPatch(i, ret);
      (*testout) << i << ": " << ret << endl;
    }

    const BitArray & freedofs_test = *fes_test->GetFreeDofs();
    ParallelForRange (Range(npatch), [&] (IntRange r)
    {
      Array<DofId> patchdofs_trial, patchdofs_test,
                    patchdofs_trial2, patchdofs_test2,
                    dofs_trial, dofs_test,
                    dofs_trial2, dofs_test2,
                    el2patch_test, el2patch_trial;
      Array<size_t> els_in_patch, facets_in_patch;
      LocalHeap lh = clh.Split ();
      for (auto p : r)
      {
        elagg->GetPatch(p, els_in_patch);
        elagg->GetInnerPatchFacets(p, facets_in_patch);
        
        patchdofs_test.SetSize0();
        patchdofs_trial.SetSize0();
        for (size_t el : els_in_patch)
        {
          Array<DofId> dofs;
          fes_trial->GetDofNrs(ElementId(VOL, el), dofs);
          for (auto d : dofs)
            if (freedofs_trial.Test(d) && !patchdofs_trial.Contains(d))
            {
              patchdofs_trial.Append(d);
            }
          fes_test->GetDofNrs(ElementId(VOL, el), dofs);
          for (auto d : dofs)
            if (freedofs_test.Test(d) && !patchdofs_test.Contains(d))
            {
              patchdofs_test.Append(d);
            }
        }

        Array<shared_ptr<BilinearFormIntegrator>> bfis[4]; // VOL, BND, ..
        Array<shared_ptr<BilinearFormIntegrator>> bfis_skeleton[4];
        for (auto icf : bf->icfs)
        {
          shared_ptr<BilinearFormIntegrator> bfi = icf->MakeBilinearFormIntegrator();
          auto fbfi = dynamic_pointer_cast<FacetBilinearFormIntegrator> (bfi);
          auto & dx = icf->dx;
          if (!fbfi)
            bfis[dx.vb] += icf->MakeBilinearFormIntegrator ();
          else
            bfis_skeleton[dx.vb] += icf->MakeBilinearFormIntegrator();
        }        

        Array<shared_ptr<LinearFormIntegrator>> lfis[4];
        if (lf)
          for (auto icf : lf->icfs)
          {
            auto &dx = icf->dx;
            lfis[dx.vb] += icf->MakeLinearFormIntegrator ();
          }


        FlatMatrix<> patchmat(patchdofs_test.Size(), patchdofs_trial.Size(), lh);
        FlatVector<> patchvec(patchdofs_test.Size(), lh);        
        patchmat = 0.0;
        patchvec = 0.0;

        for (size_t el : els_in_patch)
        {
          HeapReset hr(lh);
          ElementId ei(VOL, el);
          fes_trial->GetDofNrs(ei, dofs_trial);
          fes_test->GetDofNrs(ei, dofs_test);
          auto & trafo = ma->GetTrafo(ei, lh);
          auto & fel_trial = fes_trial->GetFE(ei, lh);
          auto & fel_test = fes_test->GetFE(ei, lh);

          MixedFiniteElement fel(fel_trial, fel_test);

          el2patch_trial.SetSize(dofs_trial.Size());
          for (auto i : Range(dofs_trial))
            el2patch_trial[i] = patchdofs_trial.Pos(dofs_trial[i]);          
          el2patch_test.SetSize(dofs_test.Size());
          for (auto i : Range(dofs_test))
            el2patch_test[i] = patchdofs_test.Pos(dofs_test[i]);   

          FlatMatrix<> elmat(dofs_trial.Size(),dofs_test.Size(), lh);
          FlatMatrix<> elmati(dofs_trial.Size(),dofs_test.Size(), lh);

          elmat = 0.0;
          for (auto & bfi : bfis[VOL])
          {
            auto & mapped_trafo = trafo.AddDeformation(bfi->GetDeformation().get(), lh);
            bfi -> CalcElementMatrix(fel, mapped_trafo, elmati, lh);
            elmat += elmati;
            // bfi -> CalcElementMatrixAdd(fel, trafo, elmat, lh);
          }        

          for (auto i : Range(dofs_test))
            for (auto j : Range(dofs_trial))
              if (el2patch_test[i] != Array<DofId>::ILLEGAL_POSITION &&
                  el2patch_trial[j] != Array<DofId>::ILLEGAL_POSITION)
                patchmat(el2patch_test[i], el2patch_trial[j]) += elmat(i,j);
          if (lf)
          {
            FlatVector<> elvec(dofs_test.Size(), lh);
            FlatVector<> sumelvec(dofs_test.Size(), lh);
            sumelvec = 0.0;
            for (auto &lfi : lfis[VOL])
            {
              auto & mapped_trafo = trafo.AddDeformation(lfi->GetDeformation().get(), lh);
              lfi->CalcElementVector(fel_test, mapped_trafo, elvec, lh);
              sumelvec += elvec;
            }

            for (auto i : Range(dofs_test))
              if (el2patch_test[i] != Array<DofId>::ILLEGAL_POSITION)
                patchvec(el2patch_test[i]) += sumelvec(i);
          }
        }
        
        for (size_t fac : facets_in_patch)
        {
          HeapReset hr(lh);
          Array<int> fac_els;
          ma->GetFacetElements(fac, fac_els);
          
          ElementId ei(VOL, fac_els[0]);
          ElementId ei2(VOL, fac_els[1]);
          fes_trial->GetDofNrs(ei, dofs_trial);
          fes_test->GetDofNrs(ei, dofs_test);
          fes_trial->GetDofNrs(ei2, dofs_trial2);
          fes_test->GetDofNrs(ei2, dofs_test2);
          auto & trafo = ma->GetTrafo(ei, lh);
          auto & trafo2 = ma->GetTrafo(ei2, lh);
          auto & fel_trial = fes_trial->GetFE(ei, lh);
          auto & fel_test = fes_test->GetFE(ei, lh);
          auto & fel_trial2 = fes_trial->GetFE(ei2, lh);
          auto & fel_test2 = fes_test->GetFE(ei2, lh);

          MixedFiniteElement fel(fel_trial, fel_test);
          MixedFiniteElement fel2(fel_trial2, fel_test2);

          el2patch_trial.SetSize(dofs_trial.Size()+dofs_trial2.Size());
          for (auto i : Range(dofs_trial.Size()+dofs_trial2.Size()) )
          {
            DofId d = i<dofs_trial.Size() ? dofs_trial[i] : dofs_trial2[i-dofs_trial.Size()];
            el2patch_trial[i] = patchdofs_trial.Pos(d); 
          }     
          el2patch_test.SetSize(dofs_test.Size()+dofs_test2.Size());
          for (auto i : Range(dofs_test.Size()+dofs_test2.Size()) )
          {
            DofId d = i<dofs_test.Size() ? dofs_test[i] : dofs_test2[i-dofs_test.Size()];
            el2patch_test[i] = patchdofs_test.Pos(d); 
          }     

          FlatMatrix<> elmat(dofs_trial.Size()+dofs_trial2.Size(),dofs_test.Size()+dofs_test2.Size(), lh);
          FlatMatrix<> elmati(dofs_trial.Size()+dofs_trial2.Size(),dofs_test.Size()+dofs_test2.Size(), lh);

          auto fnums = ma->GetElFacets(ei);
          int facnr1 = fnums.Pos(fac);             
          fnums = ma->GetElFacets(ei2);
          int facnr2 = fnums.Pos(fac);

          Array<int> vnums; vnums = ma->GetElVertices (ei);
          Array<int> vnums2; vnums2 = ma->GetElVertices (ei2);

          elmat = 0.0;
          for (auto & bfi : bfis_skeleton[VOL])
          {
            auto fbfi = dynamic_pointer_cast<FacetBilinearFormIntegrator>(bfi);
            if (fbfi) {
              auto & mapped_trafo = trafo.AddDeformation(bfi->GetDeformation().get(), lh);
              auto & mapped_trafo2 = trafo2.AddDeformation(bfi->GetDeformation().get(), lh);
              fbfi -> CalcFacetMatrix(fel, facnr1, mapped_trafo, vnums,
                                      fel2, facnr2, mapped_trafo2, vnums2, elmati, lh);
              elmat += elmati;
            }
            else throw Exception("Cast failed");
            // bfi -> CalcElementMatrixAdd(fel, trafo, elmat, lh);
          }  

          for (auto i : Range(dofs_test.Size() + dofs_test2.Size()))
            for (auto j : Range(dofs_trial.Size() + dofs_trial2.Size()))
              if (el2patch_test[i] != Array<DofId>::ILLEGAL_POSITION &&
                  el2patch_trial[j] != Array<DofId>::ILLEGAL_POSITION)
                patchmat(el2patch_test[i], el2patch_trial[j]) += elmat(i,j);
        }

        (*testout) << "patchmat in Patchloop\n" << patchmat << endl;

        patchwise_func(p,patchmat,patchvec,patchdofs_trial,patchdofs_test,lh);

      }
    });
  }

  void PatchwiseSolve(shared_ptr<ElementAggregation> elagg, 
                      shared_ptr<FESpace> fes,
                      shared_ptr<SumOfIntegrals> bf,
                      shared_ptr<SumOfIntegrals> lf,
                      shared_ptr<BaseVector> vec,
                      LocalHeap & clh
                      )
  {
    const BitArray & freedofs = *fes->GetFreeDofs();
    shared_ptr<MeshAccess> ma = elagg->GetMesh();

    size_t n_patch = elagg->GetNPatches(true);
    Array<int> dof_in_npatches(fes->GetNDof());
    dof_in_npatches = 0;
    for (int pi : Range(n_patch))
    {
      Array<DofId> patchdofs;
      Array<size_t> els_in_patch;
      patchdofs.SetSize0();
      elagg->GetPatch(pi, els_in_patch);
      //cout << pi << " - els_in_patch - " << els_in_patch << endl;
      for (int i : els_in_patch)
      {
        Array<DofId> dofs;
        fes->GetDofNrs(ElementId(VOL, i), dofs);
        for (auto d : dofs)
          if (!patchdofs.Contains(d))
          {
            patchdofs.Append(d);
            dof_in_npatches[d] ++;
          }
      }
    }    
    //cout << " - dof_in_npatches - " << dof_in_npatches << endl;
    (*vec) = 0.0;
    PatchLoop(elagg, true, fes, fes, bf, lf, clh, [&] (int p, 
                                                       FlatMatrix<> patchmat, 
                                                       FlatVector<> patchvec,
                                                       Array<DofId> & patchdofs,
                                                       Array<DofId> & patchdofs_dummy,
                                                       LocalHeap & lh
                                                       )
    {
      FlatMatrix<> patchinv(patchdofs.Size(),patchdofs.Size(),lh); 
      CalcInverse(patchmat, patchinv);
      FlatVector<> patchsol(patchdofs.Size(),lh); 
      patchsol = patchinv * patchvec;
      vec->AddIndirect (patchdofs, patchsol);
    });

    ParallelForRange (dof_in_npatches.Size(), [&] (IntRange r)
    {
      //VectorMem<10,SCAL> fluxi(dim);
      VectorMem<10,double> fluxi(1);
      ArrayMem<int,1> dnums(1);
      for (auto i : r)
        if (dof_in_npatches[i] > 1)
        {
          dnums[0] = i;
          vec->GetIndirect (dnums, fluxi);
          fluxi /= double (dof_in_npatches[i]);
          vec->SetIndirect (dnums, fluxi);
        }
    });
    
  }

  shared_ptr<SparseMatrix<double>> SetupExtensionEmbedding(shared_ptr<ElementAggregation> elagg,
                                                           shared_ptr<FESpace> fes,
                                                           shared_ptr<SumOfIntegrals> bf,
                                                           LocalHeap &clh)
  {

    // Setup of Sparse-Matrix for Embedding
    const BitArray & freedofs = *fes->GetFreeDofs();
    shared_ptr<MeshAccess> ma = elagg->GetMesh();

    BitArray non_trivial_dofs(fes->GetNDof());
    non_trivial_dofs.Clear();
    BitArray trivial_dofs(fes->GetNDof());
    trivial_dofs.Clear();
    BitArray root_dofs(fes->GetNDof());
    root_dofs.Clear();

    Array<int> dof_in_npatches(fes->GetNDof());
    dof_in_npatches = 0;

    (*testout) << "dof_in_npatches:\n" << dof_in_npatches << endl;

    ma->IterateElements(VOL, clh, [&](auto ei, LocalHeap &mlh)
    {
      Array<DofId> dofs;
      if ( (*elagg->GetElsInNontrivialPatch())[ei.Nr()] )
      {
        const bool root = (*elagg->GetRootElements())[ei.Nr()];
        fes->GetDofNrs(ei, dofs);
        for( auto i : dofs )
        {
          non_trivial_dofs.SetBitAtomic(i);
          if (root)
            root_dofs.SetBitAtomic(i);
        }
      }
      if ( (*elagg->GetElsInTrivialPatch())[ei.Nr()] )
      {
        fes->GetDofNrs(ei, dofs);
        for( auto i : dofs )
          trivial_dofs.SetBitAtomic(i);
      }

    });

    BitArray rdofs(fes->GetNDof());
    rdofs.Clear();
    rdofs |= trivial_dofs;
    rdofs |= root_dofs;
    const int nrdofs = rdofs.NumSet(); 
    const size_t n_root_dof = root_dofs.NumSet();
    (*testout) << "n_root_dof: " << n_root_dof << endl;

    size_t n_nontriv_patch = elagg->GetNPatches(false);

    for (int pi : Range(n_nontriv_patch))
    {
      Array<DofId> patchdofs;
      Array<size_t> els_in_patch;
      patchdofs.SetSize0();
      elagg->GetPatch(pi, els_in_patch);
      for (int i : els_in_patch)
      {
        Array<DofId> dofs;
        fes->GetDofNrs(ElementId(VOL, i), dofs);
        for (auto d : dofs)
          if (!patchdofs.Contains(d) && non_trivial_dofs[d] && !root_dofs[d])
          {
            patchdofs.Append(d);
            dof_in_npatches[d] ++;
          }
      }
    }

    Array<DofId> dof2rdof(fes->GetNDof());
    dof2rdof = -1;
    Array<DofId> rdof2dof(nrdofs);
    rdof2dof = -1;

    int cnt = 0;
    for (int i : Range(fes->GetNDof()))
    {
        if (!non_trivial_dofs.Test(i) || root_dofs.Test(i))
        {
        dof2rdof[i] = cnt;
        rdof2dof[cnt] = i;
        cnt++;
        }
    }

    (*testout) << "dof2rdof \n"
                << dof2rdof << endl;
    (*testout) << "rdof2dof \n"
                << rdof2dof << endl;

    size_t npatch = elagg->GetNPatches(true);
    Table<int> table, rtable;
    TableCreator<int> creator(npatch);
    TableCreator<int> creator2(npatch);
    {
      Array<size_t> els_in_patch;
      Array<DofId> patchdofs;
      for (; !creator.Done(); creator++, creator2++)
      {
        for (int pi : Range(npatch))
        {
        patchdofs.SetSize0();
        elagg->GetPatch(pi, els_in_patch);
        for (int i : els_in_patch)
        {
          Array<DofId> dofs;
          auto ei = ElementId(VOL, i);
          fes->GetDofNrs(ei, dofs);
          for (auto d : dofs)
            if (!patchdofs.Contains(d))
              patchdofs.Append(d);
        }
        for (auto d : patchdofs)
        {
          creator.Add(pi,d);
          const int rdof = dof2rdof[d];
          if (IsRegularDof(rdof))
            creator2.Add(pi,rdof);
        }
        }
      }
      table = creator.MoveTable ();
      rtable = creator2.MoveTable ();
    }
    shared_ptr<SparseMatrix<double>> P
      = make_shared<SparseMatrix<double>> (fes->GetNDof (), 
                                           nrdofs, table, rtable, false);
    
    P->SetZero ();

    double one = 1;
    FlatMatrix<> I(1,1,&one);
    Array<DofId> cdof(1), crdof(1);
    for (int i : Range(fes->GetNDof()))
    {
      if (!non_trivial_dofs[i] || root_dofs[i])
      {
        cdof[0] = i; 
        crdof[0] = dof2rdof[i];
        P->AddElementMatrix(cdof, crdof, I);
      }
    }

    (*testout) << "non_trivial_dofs: " << non_trivial_dofs << endl;
    (*testout) << "P: \n" << *P << endl;

    PatchLoop(elagg, false, fes, fes, bf, nullptr, clh, [&] (int p, 
                                                            FlatMatrix<> patchmat, 
                                                            FlatVector<> patchvec,
                                                            Array<DofId> & patchdofs,
                                                            Array<DofId> & patchdofs_dummy,
                                                            LocalHeap & lh
                                                            )
    {
      Array<size_t> els_in_patch;
      elagg->GetPatch(p, els_in_patch);

      Array<DofId> dofs;
      fes->GetDofNrs(ElementId(VOL, els_in_patch[0]), dofs);
      int nr_dofs_inner = dofs.Size();

      /*int nr_dofs_inner = 0;
      for (auto d : dofs)
        if (freedofs.Test(d))
            nr_dofs_inner++;*/

      (*testout) << "nr dofs root " << nr_dofs_inner << endl;

      int nr_dofs_outer = patchdofs.Size() - nr_dofs_inner;

      FlatMatrix<double> Soo(nr_dofs_outer, nr_dofs_outer, lh);
      Soo = patchmat.Rows(nr_dofs_inner, patchdofs.Size()).Cols(nr_dofs_inner, patchdofs.Size());

      FlatMatrix<double> Soi(nr_dofs_outer, nr_dofs_inner, lh);
      Soi = patchmat.Rows(nr_dofs_inner, patchdofs.Size()).Cols(0, nr_dofs_inner);

      FlatMatrix<double> Sooinv(nr_dofs_outer, nr_dofs_outer, lh);
      CalcInverse(Soo, Sooinv);
      FlatMatrix<double> Q(nr_dofs_outer, nr_dofs_inner, lh);
      Q = -Sooinv * Soi;

      for (int i : Range(nr_dofs_outer))
      {
        Q.Row(i) /= dof_in_npatches[patchdofs[nr_dofs_inner+i]];
      }
      Array<DofId> patchrdofs(nr_dofs_inner);
      for (int i : Range(nr_dofs_inner))
        patchrdofs[i] = dof2rdof[patchdofs[i]];
      (*testout) << "p = " << p << "\n patchdofs: \n" << patchdofs << "\ninner dofs: \n" << patchrdofs << endl;
      P->AddElementMatrix(patchdofs.Range(nr_dofs_inner,patchdofs.Size()), patchrdofs, Q);

      (*testout) << "patch nr " << p << endl;
      // (*testout) << "dofs test\n" << patchdofs_test << endl;
      // (*testout) << "patch vector \n" << patchvec << endl;
      (*testout) << "patch matrix \n"
                 << patchmat << endl;
      (*testout) << "patch dofs\n"
                 << patchdofs << endl;
      // (*testout) << "dofs of root \n" << dofs << endl;
      (*testout) << "Q: \n"
                 << Q << endl;

    });

    // Finalizing steps
    (*testout) << "P result: \n"
               << *P << endl;  
    return P;
  } // end of SetupAgg..
}