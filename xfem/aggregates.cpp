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

  }


  // This is a dummy as a first step towards 
  //    * solutions of patchwise problems 
  //    * or setup of embedding matrices
  //    * or similar things
  void PatchDummy (shared_ptr<ElementAggregation> elagg, 
                   shared_ptr<FESpace> fes_trial,
                   shared_ptr<FESpace> fes_test,
                   shared_ptr<SumOfIntegrals> bf,
                   shared_ptr<SumOfIntegrals> lf
                  )
  {
    cout << " hello from dummy " << endl;
    //size_t npatch = 

  }



}
