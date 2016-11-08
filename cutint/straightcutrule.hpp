#pragma once
#include "xintegration.hpp"
#include <memory>
#include <numeric>

using namespace ngfem;
namespace xintegration
{
  typedef Array<tuple<Vec<3>, double>> PointCnt;
  DOMAIN_TYPE CheckIfStraightCut(FlatVector<> cf_lset_at_element);

  class Polytope {
  public:
      Array<int> ia; //the points
      int D; //Dimension
      shared_ptr<PointCnt> svs_ptr;
      Polytope(Array<int> & a_ia, int a_D, auto a_svs_ptr) : ia(a_ia), D(a_D), svs_ptr(a_svs_ptr) {;}
      Polytope(initializer_list<int> a_ia_l, int a_D, auto a_svs_ptr) : ia(a_ia_l), D(a_D), svs_ptr(a_svs_ptr) {;}
      Polytope(){D = -1; svs_ptr = nullptr; }

      auto begin() {return ia.begin(); }
      auto end() {return ia.end();}
      auto begin() const {return ia.begin(); }
      auto end() const {return ia.end();}

      void Append(auto v){ ia.Append(v); }

      auto operator[](int j) const { return ia[j]; }
      auto Size() const { return ia.Size(); }
      void DeleteElement(auto i){ ia.DeleteElement(i);}

      Vec<3> GetPoint(int j) const {
          return get<0>((*svs_ptr)[ia[j]]);
      }
  };

  class StraightCutElementGeometry {      
  private:
      int D;
      shared_ptr<PointCnt> svs_ptr;
      Array<Polytope> simplices;
      FlatVector<> lset;
      ELEMENT_TYPE et;
      double MeasureSimplVol(const Polytope &s);
      Polytope CalcCutPolytopeUsingLset(const Polytope &s);
      Polytope CalcCutPointLineUsingLset(const Polytope &s);
      void CalcNormal();
      void LoadBaseSimplexFromElementTopology();
      void CutBaseSimplex(DOMAIN_TYPE dt);

  public:
      LocalHeap & lh;
      Vec<3> normal;

      StraightCutElementGeometry(FlatVector<> a_lset, ELEMENT_TYPE a_et, LocalHeap &a_lh) : lset(a_lset), et(a_et), lh(a_lh) {
          D = Dim(et);
          svs_ptr = make_shared<PointCnt>();
      }

      void GetIntegrationRule(int order, DOMAIN_TYPE dt, IntegrationRule &intrule);
  };

  class StraightCutQuadElementGeometry {
  public:
      int D;
      ELEMENT_TYPE et;
      Polytope base_quad;
  };

  const IntegrationRule * StraightCutIntegrationRule(shared_ptr<CoefficientFunction> cf_lset,
                                                     const FlatVector<> & cf_lset_at_element,
                                                     const ElementTransformation & trafo,
                                                     DOMAIN_TYPE dt,
                                                     int intorder,
                                                     LocalHeap & lh);
}
