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
    std::cout << "Hello from element aggregation" << endl;
  }

  void ElementAggregation::Update(shared_ptr<BitArray> & root_els, shared_ptr<BitArray> & bad_els, LocalHeap & lh){
    const int ne = ma->GetNE();

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
    shared_ptr<BitArray>cluster_root_elements = make_shared<BitArray>(ne);
    cluster_root_elements->Clear();
    //cout << *cluster_root_elements << endl;
    *cluster_root_elements |= *root_els;
    //cout << *cluster_root_elements << endl;
    *cluster_root_elements &= *GetElementsWithNeighborFacets(ma,facets,lh);
    //cout << *cluster_root_elements << endl;
    int nc = 0;
    for (int i : Range(ne)){
      if (cluster_root_elements->Test(i)){
        nc ++;
      }
    }
    cout << "Number of cluster: " << nc << endl;

    Vector<> element_to_cluster(ne);
    element_to_cluster = 0.0;

    Array<int>cluster_roots(nc);
    int counter = 0;
    for (int i: Range(ma->GetNE())){
      if (cluster_root_elements->Test(i)){
        element_to_cluster(i) = counter;
        cluster_roots[counter]=i;
        counter++;
      }
      else
        element_to_cluster(i) = -1.0;
    }

    BitArray front_elements(ne);
    BitArray new_front_elements(ne);
    new_front_elements = *cluster_root_elements;

    BitArray not_covered_yet(ne);
    not_covered_yet.Clear();
    not_covered_yet |= *bad_els;
    BitArray newly_covered(ne);
    bool uncovered_elements_left = true;

    while (not_covered_yet.NumSet() > 0)
    {
      front_elements = new_front_elements;
      new_front_elements.Clear();
      newly_covered.Clear();
      //IterateRange(ne, lh, [&] (int elnr, LocalHeap & lh)
      for (int elnr = 0; elnr < ne; elnr++)
      {
        int cluster_id = element_to_cluster[elnr];
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
              if (not_covered_yet.Test(nelnr) && (!newly_covered.Test(nelnr)))
              {
                new_front_elements.SetBitAtomic(nelnr);
                element_to_cluster[nelnr] = cluster_id;
                newly_covered.SetBitAtomic(nelnr);
              }
            }
          }

        }
      }
      not_covered_yet &= ~newly_covered;


    }
    cout << "element_to_cluster = \n" << element_to_cluster << endl;

 }


}
