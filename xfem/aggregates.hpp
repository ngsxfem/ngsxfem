#pragma once

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>


using namespace ngsolve;
// using namespace cutinfo;

namespace ngcomp
{


  class ElementAggregation
  {
    
  protected:
    shared_ptr<MeshAccess> ma;
    shared_ptr<BitArray> inner_patch_facets;
    Array<size_t> cluster_roots;
    Table<size_t> cluster_leafs;
    Table<size_t> cluster_facets;
  public:
    ElementAggregation (shared_ptr<MeshAccess> ama);
    
    shared_ptr<MeshAccess> GetMesh () const { return ma; }


    void Update(shared_ptr<BitArray> & root_els, shared_ptr<BitArray> & bad_els, LocalHeap & lh);
  
    shared_ptr<BitArray> & GetInnerPatchFacets(){ return inner_patch_facets; }
    
  
  };



}