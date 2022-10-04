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
  shared_ptr<BitArray>initial_facets = GetFacetsWithNeighborTypes(ma,root_els,bad_els,false,false,true,lh);
  shared_ptr<BitArray>cluster_root_elements = make_shared<BitArray>(ma->GetNE());
  cluster_root_elements->Clear();
  //cout << *cluster_root_elements << endl;
  *cluster_root_elements |= *root_els;
  //cout << *cluster_root_elements << endl;
  *cluster_root_elements &= *GetElementsWithNeighborFacets(ma,initial_facets,lh);
  //cout << *cluster_root_elements << endl;
  int nc = 0;
  for (int i : Range(ma->GetNE())){
    if (cluster_root_elements->Test(i)){
      nc ++;
    }
  }
  cout << "Number of cluster: " << nc << endl;
  Array<int>cluster_roots(nc);
  int counter = 0;
  for (int i: Range(ma->GetNE())){
    if (cluster_root_elements->Test(i)){
      cluster_roots[counter]=i;
      counter++;
    }
  }
 }
 /* 
need: array front; array front_to_cluster; BitArray not_covered_yet
need: cluster_table with suitable type / size
loop until not_covered_yet is empty:
  front = new_front
  front_to_cluster = new_front_to_cluster
  new_front, new_front_to_cluster = empty
  loop (element_index, cluster_index) over front and front_to_cluster:
    els = neighbors[element_index]
    loop el in els: 
      if el == Bad and el in not_covered_yet: 
        new_front.add(el) 
        new_front_to_cluster = cluster_index
        not_covered_yet(el) = False
        add element to cluster_table
 */
}
