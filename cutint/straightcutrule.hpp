#pragma once
#include "xintegration.hpp"
#include <algorithm>
#include <numeric>

using namespace ngfem;
namespace xintegration
{
  DOMAIN_TYPE CheckIfStraightCut(FlatVector<> cf_lset_at_element);

  typedef Array<int> SimpleX;

  class StraightCutElementGeometry {      
  private:
      int D;
      double MeasureSimplVol(SimpleX &s);
      SimpleX Cut(SimpleX &s);
      void AppendIntegrationRuleOnSimpl(SimpleX &s, int order, IntegrationRule &intrule);
      void CalcNormal(SimpleX &s_cut);
  public:
      vector<Vec<3>> svs; //Simplex Vertices
      vector<SimpleX> simplices;
      FlatVector<> lset;
      ELEMENT_TYPE et;
      LocalHeap & lh;
      Vec<3> normal;

      StraightCutElementGeometry(FlatVector<> a_lset, ELEMENT_TYPE a_et, LocalHeap &a_lh) : lset(a_lset), et(a_et), lh(a_lh) {
          D = Dim(et);
      }

      void LoadBaseSimplexFromElementTopology();
      void CutBaseSimplex(DOMAIN_TYPE dt);
      void GetIntegrationRule(int order, IntegrationRule &intrule);
  };

  const IntegrationRule * StraightCutIntegrationRule(shared_ptr<CoefficientFunction> cf_lset,
                                                     const FlatVector<> & cf_lset_at_element,
                                                     const ElementTransformation & trafo,
                                                     DOMAIN_TYPE dt,
                                                     int intorder,
                                                     LocalHeap & lh);
  const IntegrationRule * StraightCutIntegrationRuleOld(shared_ptr<CoefficientFunction> cf_lset,
                                                       const FlatVector<> & cf_lset_at_element,
                                                       const ElementTransformation & trafo,
                                                       DOMAIN_TYPE dt,
                                                       int intorder,
                                                       LocalHeap & lh);
}
