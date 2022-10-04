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
  int counter = 0;
  for (int i : Range(ma->GetNE())){
    if (cluster_root_elements->Test(i)){
      counter ++;
    }
  }
  cout << "Counter=" << counter << endl;
 }

}
