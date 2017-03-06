#pragma once

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

/// from ngxfem
#include "../cutint/xintegration.hpp"
// #include "../xfem/xfiniteelement.hpp"

using namespace ngsolve;
using namespace xintegration;
// using namespace cutinfo;

namespace ngcomp
{

  class XFESpace;

  class CutInformation
  {
    friend class XFESpace;
  protected:
    shared_ptr<MeshAccess> ma;
    Array<shared_ptr<VVector<double>>> cut_ratio_of_element;
    Array<shared_ptr<BitArray>> elems_of_domain_type;
    Array<shared_ptr<BitArray>> selems_of_domain_type;
    Array<shared_ptr<BitArray>> facets_of_domain_type;
    Array<shared_ptr<BitArray>> cut_neighboring_node;

    double subdivlvl = 0;
  public:
    CutInformation (shared_ptr<MeshAccess> ama);
    void Update(shared_ptr<CoefficientFunction> lset, LocalHeap & lh);

    shared_ptr<MeshAccess> GetMesh () const { return ma; }

    shared_ptr<BaseVector> GetCutRatios (VorB vb) const
    {
      return cut_ratio_of_element[vb];
    }

    // template <NODE_TYPE NT>
    // DOMAIN_TYPE GetDomainOfNode (int nr)
    // {
    //   if (GetCutRatioOfNode<NT>(nr) > 0.0 && GetCutRatioOfNode<NT>(nr) < 1.0)
    //     return IF;
    //   else
    //     return GetCutRatioOfNode<NT>(nr) == 0.0 ? NEG : POS;
    // }

    shared_ptr<BitArray> GetElementsOfDomainType(DOMAIN_TYPE dt, VorB vb) const
    {
      if (vb == VOL)
        return elems_of_domain_type[dt];
      else
        return selems_of_domain_type[dt];
    }

    shared_ptr<BitArray> GetFacetsOfDomainType(DOMAIN_TYPE dt) const { return facets_of_domain_type[dt]; }

  };

  shared_ptr<BitArray> GetFacetsWithNeighborTypes(shared_ptr<MeshAccess> ma,
                                                  shared_ptr<BitArray> a,
                                                  shared_ptr<BitArray> b,
                                                  bool bound_val_a,
                                                  bool bound_val_b,
                                                  bool ask_and,
                                                  LocalHeap & lh);
  shared_ptr<BitArray> GetElementsWithNeighborFacets(shared_ptr<MeshAccess> ma,
                                                     shared_ptr<BitArray> a,
                                                     LocalHeap & lh);

}
