#pragma once

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>
// #include <integratorcf.hpp>
#include "../cutint/cutintegral.hpp"

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
    shared_ptr<BitArray> el_roots = nullptr;
    size_t n_trivial_els = -1;
    size_t n_els_in_nontrivial_patches = -1;
    size_t n_nontrivial_patches = -1;
    Vector<int> element_to_patch;
    Vector<int> facet_to_patch;
    Array<size_t> patch_roots;
    Table<size_t> patch_leafs;
    Table<size_t> patch_facets;
    Table<size_t> patch;
    Array<int> trivial_patch_to_element;

  public:
    ElementAggregation(shared_ptr<MeshAccess> ama);

    shared_ptr<MeshAccess> GetMesh() const { return ma; }

    void Update(shared_ptr<BitArray> &root_els, shared_ptr<BitArray> &bad_els, LocalHeap &lh);

    shared_ptr<BitArray> &GetInnerPatchFacets() { return inner_patch_facets; }
    void GetInnerPatchFacets(int patchnr, Array<size_t> &ret);
    shared_ptr<BitArray> &GetElsInTrivialPatch() { return el_in_trivial_patch; }
    shared_ptr<BitArray> &GetElsInNontrivialPatch() { return el_in_nontrivial_patch; }
    shared_ptr<BitArray> &GetRootElements() { return el_roots; }
    const Array<size_t> &GetPatchRoots() { return patch_roots; }
    const Table<size_t> &GetPatchLeafs() { return patch_leafs; }
    const Table<size_t> &GetPatchFacets() { return patch_facets; }
    const Table<size_t> &GetPatch() { return patch; }
    void GetPatch(int patchnr, Array<size_t> &ret);
    const Vector<int> &GetElementToPatch() { return element_to_patch; }
    const Vector<int> &GetFacetToPatch() { return facet_to_patch; }

    size_t GetNNontrivialPatches() { return n_nontrivial_patches; }
    size_t GetNPatches(bool count_trivial) { return n_nontrivial_patches + (count_trivial ? n_trivial_els : 0); }
  };



  template <typename TFUNC>
  void PatchLoop (shared_ptr<ElementAggregation> elagg, 
                  bool include_trivial_patches,
                  shared_ptr<FESpace> fes_trial,
                  shared_ptr<FESpace> fes_test,
                  shared_ptr<SumOfIntegrals> bf,
                  shared_ptr<SumOfIntegrals> lf,
                  LocalHeap & clh,
                  const TFUNC & patchwise_func
                  );

  // setup AggFEM Embedding matrix
  shared_ptr<SparseMatrix<double>> SetupExtensionEmbedding(shared_ptr<ElementAggregation> elagg,
                                                     shared_ptr<FESpace> fes,
                                                     shared_ptr<SumOfIntegrals> bf,
                                                     LocalHeap &clh);

  // solve patch-local problems
  void PatchwiseSolve(shared_ptr<ElementAggregation> elagg, 
                      shared_ptr<FESpace> fes,
                      shared_ptr<SumOfIntegrals> bf,
                      shared_ptr<SumOfIntegrals> lf,
                      shared_ptr<BaseVector> vec,
                      LocalHeap & clh
                      );
}
