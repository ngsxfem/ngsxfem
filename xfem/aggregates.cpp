/// from ngxfem
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

  void ElementAggregation::Update(shared_ptr<BitArray> & root_els, shared_ptr<BitArray> & bad_els, LocalHeap & lh){
    const size_t ne = ma->GetNE();
    const size_t nf = ma->GetNFacets();

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
    cout << "Number of (nontrivial) patches: " << nc << endl;

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

  void ElementAggregation::GetPatch(int patchnr, Array<size_t> & ret){
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

  void ElementAggregation::GetInnerPatchFacets(int patchnr, Array<size_t> & ret){
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



  // This is a dummy as a first step towards 
  //    * solutions of patchwise problems ≈
  //    * or setup of embedding matrices
  //    * or similar things
  void PatchDummy (shared_ptr<ElementAggregation> elagg, 
                   shared_ptr<FESpace> fes_trial,
                   shared_ptr<FESpace> fes_test,
                   shared_ptr<SumOfIntegrals> bf,
                   shared_ptr<SumOfIntegrals> lf,
                   LocalHeap & clh
                  )
  {
    cout << " hello from dummy " << endl;
    shared_ptr<MeshAccess> ma = elagg->GetMesh();
    Array<size_t> ret;
    for (int i : Range(elagg->GetNPatches( /*count_trivial=*/false )) )
    {
      elagg->GetPatch(i, ret);
      cout << i << ": " << ret << endl;
    }

    size_t npatch = elagg->GetNPatches(/*count_trivial=*/false);

    const BitArray & freedofs_test = *fes_test->GetFreeDofs();
    const BitArray & freedofs_trial = *fes_trial->GetFreeDofs();
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
                patchdofs_trial.Append(d);
          fes_test->GetDofNrs(ElementId(VOL, el), dofs);
          for (auto d : dofs)
            if (freedofs_test.Test(d) && !patchdofs_test.Contains(d))
                patchdofs_test.Append(d);
        }

        Array<shared_ptr<BilinearFormIntegrator>> bfis[4]; // VOL, BND, ..
        Array<shared_ptr<BilinearFormIntegrator>> bfis_skeleton[4];
        for (auto icf : bf->icfs)
        {
          auto & dx = icf->dx;
          if(!dx.skeleton && 0!=dynamic_pointer_cast<FacetPatchDifferentialSymbol>(make_shared<DifferentialSymbol>(dx)))
            bfis[dx.vb] += icf->MakeBilinearFormIntegrator ();
          else
            bfis_skeleton[dx.vb] += icf->MakeBilinearFormIntegrator();
          //if (dx.vb == VOL)
          //  bfis += make_shared<SymbolicBilinearFormIntegrator> (icf->cf, VOL, dx.element_vb);
        }        

        Array<shared_ptr<LinearFormIntegrator>> lfis[4];
        for (auto icf : lf->icfs)
        {
          auto &dx = icf->dx;
          lfis[dx.vb] += icf->MakeLinearFormIntegrator ();
        }


        FlatMatrix<> patchmat(patchdofs_test.Size(), patchdofs_trial.Size(), lh);
        FlatVector<> patchvec(patchdofs_test.Size(), lh);        
        //FlatVector<> patchsol(patchdofs_trial.Size(), lh);        
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
            bfi -> CalcElementMatrix(fel, trafo, elmati, lh);
            elmat += elmati;
            // bfi -> CalcElementMatrixAdd(fel, trafo, elmat, lh);
          }        

            
          FlatVector<> elvec(dofs_test.Size(), lh);
          FlatVector<> sumelvec(dofs_test.Size(), lh);
          sumelvec = 0.0;
          for (auto & lfi : lfis[VOL])
          {
            lfi -> CalcElementVector(fel_test, trafo, elvec, lh);
            sumelvec += elvec;
          }    

          for (auto i : Range(dofs_test))
            for (auto j : Range(dofs_trial))
              if (el2patch_test[i] != Array<DofId>::ILLEGAL_POSITION &&
                  el2patch_trial[j] != Array<DofId>::ILLEGAL_POSITION)
                patchmat(el2patch_test[i], el2patch_trial[j]) += elmat(i,j);
          for (auto i : Range(dofs_test))
            if (el2patch_test[i] != Array<DofId>::ILLEGAL_POSITION)
              patchvec(el2patch_test[i]) += sumelvec(i);

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
              fbfi -> CalcFacetMatrix(fel, facnr1, trafo, vnums,
                                      fel2, facnr2, trafo2, vnums2, elmati, lh);
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
        

        (*testout) << "patch nr " << p << endl;
        // (*testout) << "dofs trial \n" << patchdofs_trial << endl;
        // (*testout) << "dofs test\n" << patchdofs_test << endl;
        // (*testout) << "patch vector \n" << patchvec << endl;
        (*testout) << "patch matrix \n" << patchmat << endl;
        // TODO : + Facet loop       

      }

    } );  
  }

}
