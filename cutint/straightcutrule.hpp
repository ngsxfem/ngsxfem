#pragma once
#include "xintegration.hpp"
#include <algorithm>
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
    Polytope(Array<int> & a_ia, int a_D, shared_ptr<PointCnt> a_svs_ptr) : ia(a_ia), D(a_D), svs_ptr(a_svs_ptr) {;}
    Polytope(initializer_list<int> a_ia_l, int a_D, shared_ptr<PointCnt> a_svs_ptr) : ia(a_ia_l), D(a_D), svs_ptr(a_svs_ptr) {;}
      template<int D>
    Polytope(Vec<D, tuple<Vec<3>, double>> a_points, int a_D, shared_ptr<PointCnt> a_svs_ptr) : D(a_D), svs_ptr(a_svs_ptr)
    {
          for(int i=0; i<D; i++){
              auto p = a_points[i];
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

    void Append(int v){ ia.Append(v); }

      auto operator[](int j) const { return ia[j]; }
      auto Size() const { return ia.Size(); }
    void DeleteElement(int i){ ia.DeleteElement(i);}

      Vec<3> GetPoint(int j) const { return get<0>((*svs_ptr)[ia[j]]); }
      double GetLset(int j) const{ return get<1>((*svs_ptr)[ia[j]]); }
  };

  double MeasureSimplVol(const Polytope &s);
  Polytope CalcCutPolytopeUsingLset(const Polytope &s);
  Polytope CalcCutPointLineUsingLset(const Polytope &s);

  class CutElementGeometry {
  protected:
      ELEMENT_TYPE et;
      int D;
  public:
      virtual Vec<3> GetNormal(const Vec<3>& p) const = 0;
      virtual void GetIntegrationRule(int order, DOMAIN_TYPE dt, IntegrationRule &intrule) = 0;
  };

  class CutSimplexElementGeometry : public CutElementGeometry {
  private:
      shared_ptr<PointCnt> svs_ptr;
      Array<Polytope> simplices;
      FlatVector<> lset;

      void CalcNormal(const Polytope &base_simplex);
      void LoadBaseSimplexFromElementTopology();
      void CutBaseSimplex(DOMAIN_TYPE dt);
      Vec<3> normal;

  public:
      LocalHeap & lh;

      CutSimplexElementGeometry(FlatVector<> a_lset, ELEMENT_TYPE a_et, LocalHeap &a_lh) : lset(a_lset), lh(a_lh) {
          et = a_et;
          D = Dim(et);
          svs_ptr = make_shared<PointCnt>();
      }

      virtual Vec<3> GetNormal(const Vec<3>& p) const{
          return normal;
      }

      virtual void GetIntegrationRule(int order, DOMAIN_TYPE dt, IntegrationRule &intrule);
  };

  class CutQuadElementGeometry : public CutElementGeometry {
  private:
      shared_ptr<PointCnt> svs_ptr;
      Vec<2, Array<Polytope>> segments;
      Vec<2, Vec<2, Vec<2, double>>> lc; //TODO: Replace with MultilinearFunction
      Array<Polytope> Cut_quads; Array<Polytope> Volume_quads;
      Vector<> lset;

      void IntegrateCutQuads(int order, IntegrationRule &intrule);
      void IntegrateCutQuads3D(int order, IntegrationRule &intrule);
      void IntegrateVolumeQuads(int order, IntegrationRule &intrule);
      void IntegrateVolumeQuads3D(int order, IntegrationRule &intrule);
      void IntegrateVolumeOfCutQuads(DOMAIN_TYPE dt, int order, IntegrationRule &intrule);
      void IntegrateVolumeOfCutQuads3D(DOMAIN_TYPE dt, int order, IntegrationRule &intrule);
      void FindVolumeAndCutQuads(DOMAIN_TYPE dt);
      void FindVolumeAndCutQuads3D(DOMAIN_TYPE dt);
  public:
      void LoadBaseQuadFromElementTopology();
      LocalHeap & lh;

      virtual Vec<3> GetNormal(const Vec<3>& p) const;
      CutQuadElementGeometry(FlatVector<> a_lset, ELEMENT_TYPE a_et, LocalHeap &a_lh) : lset(a_lset), lh(a_lh) {
          et = a_et;
          D = Dim(et); svs_ptr = make_shared<PointCnt>();
      }

      virtual void GetIntegrationRule(int order, DOMAIN_TYPE dt, IntegrationRule &intrule);
  };

  const IntegrationRule * StraightCutIntegrationRule(const FlatVector<> & cf_lset_at_element,
                                                     const ElementTransformation & trafo,
                                                     DOMAIN_TYPE dt,
                                                     int intorder,
                                                     LocalHeap & lh);
}
