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
    /// @brief elements that are good (a.k.a. root), but not involved in non-trivial patch
    shared_ptr<BitArray> el_in_trivial_patch; 
    /// @brief elements that are involved in a non-trivial patches
    shared_ptr<BitArray> el_in_nontrivial_patch;
    size_t n_trivial_els = -1;
    size_t n_els_in_nontrivial_patches = -1;
    size_t n_nontrivial_patches = -1;
    Array<size_t> patch_roots;
    Table<size_t> patch_leafs;
    Table<size_t> patch_facets;
  public:
    ElementAggregation (shared_ptr<MeshAccess> ama);
    
    shared_ptr<MeshAccess> GetMesh () const { return ma; }


    void Update(shared_ptr<BitArray> & root_els, shared_ptr<BitArray> & bad_els, LocalHeap & lh);
  
    shared_ptr<BitArray> & GetInnerPatchFacets(){ return inner_patch_facets; }
    shared_ptr<BitArray> & GetElsInTrivialPatch(){ return el_in_trivial_patch; }
    shared_ptr<BitArray> & GetElsInNontrivialPatch(){ return el_in_nontrivial_patch; }
    const Array<size_t> & GetPatchRoots(){ return patch_roots; }
    const Table<size_t> & GetPatchLeafs(){ return patch_leafs; }
    const Table<size_t> & GetPatchFacets(){ return patch_facets; }
    
    size_t GetNNontrivialPatches() { return n_nontrivial_patches; }
  
  };


  // This is a dummy as a first step towards 
  //    * solutions of patchwise problems 
  //    * or setup of embedding matrices
  //    * or similar things
  void PatchDummy (shared_ptr<ElementAggregation>, 
                   shared_ptr<FESpace> fes_trial,
                   shared_ptr<FESpace> fes_test,
                   shared_ptr<SumOfIntegrals> bf,
                   shared_ptr<SumOfIntegrals> lf
                  );

}
