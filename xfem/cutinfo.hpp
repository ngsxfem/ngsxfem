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

  class CutInformation
  {
  protected:
    shared_ptr<MeshAccess> ma;
    Array<shared_ptr<VVector<double>>> cut_ratio_of_node;
    Array<shared_ptr<BitArray>> elems_of_domain_type;
    Array<shared_ptr<BitArray>> facets_of_domain_type;
  public:
    CutInformation (shared_ptr<MeshAccess> ama);
    void Update(shared_ptr<CoefficientFunction> lset);

    template <NODE_TYPE NT>
    double GetCutRatioOfNode (int nr) const
    {
      return (*cut_ratio_of_node[NT])(nr);
    }

    template <NODE_TYPE NT>
    DOMAIN_TYPE GetDomainOfNode (int nr)
    {
      if (GetCutRatioOfNode<NT>(nr) > 0.0 && GetCutRatioOfNode<NT>(nr) < 1.0)
        return IF;
      else
        return GetCutRatioOfNode<NT>(nr) == 0.0 ? NEG : POS;
    }

    shared_ptr<BitArray> GetElementsOfDomainType(DOMAIN_TYPE dt) const { return elems_of_domain_type[dt]; }
    shared_ptr<BitArray> GetFacetsOfDomainType(DOMAIN_TYPE dt) const { return facets_of_domain_type[dt]; }

  };
}
