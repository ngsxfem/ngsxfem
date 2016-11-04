#pragma once
#include "xintegration.hpp"
#include <algorithm>
#include <numeric>

using namespace ngfem;
namespace xintegration
{
  DOMAIN_TYPE CheckIfStraightCut(FlatVector<> cf_lset_at_element);

  class Polytope {
  public:
      Array<int> ia; //the points
      int D; //Dimension
      Polytope(Array<int> & a_ia, int a_D) : ia(a_ia), D(a_D) {;}
      Polytope(initializer_list<int> a_ia_l, int a_D) : ia(a_ia_l), D(a_D) {;}
      Polytope(){D = -1; }

      auto begin() {return ia.begin(); }
      auto end() {return ia.end();}
      auto begin() const {return ia.begin(); }
      auto end() const {return ia.end();}

      void Append(auto v){ ia.Append(v); }

      auto operator[](int j) const { return ia[j]; }
      auto Size() const { return ia.Size(); }
      void DeleteElement(auto i){ ia.DeleteElement(i); }
  };

  class StraightCutElementGeometry {      
  public:      int D;
      double MeasureSimplVol(const Polytope &s);
      Polytope CalcCutPolytopeUsingLset(const Polytope &s);
      Polytope CalcCutPointLineUsingLset(const Polytope &s);
      //void CalcNormal();
      void CalcNormal();
  public:
      Array<Vec<3>> svs;
      Array<Polytope> simplices;
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
