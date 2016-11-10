#pragma once
#include "xintegration.hpp"
#include <memory>
#include <numeric>

using namespace ngfem;

namespace ngbla {
bool operator==(const Vec<3> a, const Vec<3> b){
    return L2Norm(a - b)<1e-12;
}
}


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
      Polytope(initializer_list<tuple<Vec<3>, double>> a_points, int a_D, auto a_svs_ptr) : D(a_D), svs_ptr(a_svs_ptr){
          for(auto p : a_points){
              if(svs_ptr->Contains(p)){
                  ia.Append(svs_ptr->Pos(p));
              }
              else {
                  svs_ptr->Append(p);
                  ia.Append(svs_ptr->Size() - 1);
              }
          }
      }

      Polytope(){D = -1; svs_ptr = nullptr; }

      auto begin() {return ia.begin(); }
      auto end() {return ia.end();}
      auto begin() const {return ia.begin(); }
      auto end() const {return ia.end();}

      void Append(auto v){ ia.Append(v); }

      auto operator[](int j) const { return ia[j]; }
      auto Size() const { return ia.Size(); }
      void DeleteElement(auto i){ ia.DeleteElement(i);}

      Vec<3> GetPoint(int j) const { return get<0>((*svs_ptr)[ia[j]]); }
      double GetLset(int j) const{ return get<1>((*svs_ptr)[ia[j]]); }
  };

  double MeasureSimplVol(const Polytope &s);
  Polytope CalcCutPolytopeUsingLset(const Polytope &s);
  Polytope CalcCutPointLineUsingLset(const Polytope &s);

  class StraightCutElementGeometry {      
  private:
      int D;
      shared_ptr<PointCnt> svs_ptr;
      Array<Polytope> simplices;
      FlatVector<> lset;
      ELEMENT_TYPE et;
      void CalcNormal(const Polytope &base_simplex);
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
      shared_ptr<PointCnt> svs_ptr;
      Polytope base_quad;
      Array<Polytope> segments_x;
      Array<Polytope> segments_y;
      double a,b,c,d;

      void PartitionSegmentsX(double x_cut, double lset_on_x_cut);
      void PartitionSegmentsY(double y_cut, double lset_on_y_cut);

      StraightCutQuadElementGeometry(ELEMENT_TYPE a_et) : et(a_et) {D = Dim(et); svs_ptr = make_shared<PointCnt>();}

      void LoadBaseSimplexFromElementTopology(FlatVector<> lset); //TODO:Rename
      void GetIntegrationRule(FlatVector<> lset, int order, DOMAIN_TYPE dt, IntegrationRule &intrule);
  };

  const IntegrationRule * StraightCutIntegrationRule(shared_ptr<CoefficientFunction> cf_lset,
                                                     const FlatVector<> & cf_lset_at_element,
                                                     const ElementTransformation & trafo,
                                                     DOMAIN_TYPE dt,
                                                     int intorder,
                                                     LocalHeap & lh);
}
