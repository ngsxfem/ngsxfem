#pragma once
#include "xintegration.hpp"
#include <algorithm>
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
      template<int D>
      Polytope(Vec<D, tuple<Vec<3>, double>> a_points, int a_D, auto a_svs_ptr) : D(a_D), svs_ptr(a_svs_ptr){
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
      Vec<3> normal;

  public:
      LocalHeap & lh;

      Vec<3> GetNormal(const Vec<3>& p) const{
          return normal;
      }

      StraightCutElementGeometry(FlatVector<> a_lset, ELEMENT_TYPE a_et, LocalHeap &a_lh) : lset(a_lset), et(a_et), lh(a_lh) {
          D = Dim(et);
          svs_ptr = make_shared<PointCnt>();
      }

      void GetIntegrationRule(int order, DOMAIN_TYPE dt, IntegrationRule &intrule);
  };

  class StraightCutQuadElementGeometry {
  private:
      int D;
      ELEMENT_TYPE et;
      shared_ptr<PointCnt> svs_ptr;
      Vec<2, Array<Polytope>> segments;
      Vec<2, Vec<2, Vec<2, double>>> lc;
      Array<Polytope> Cut_quads; Array<Polytope> Volume_quads;
      FlatVector<> lset;

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
  public:
      Vec<3> GetNormal(const Vec<3>& p) const;
      StraightCutQuadElementGeometry(FlatVector<> a_lset, ELEMENT_TYPE a_et, LocalHeap &a_lh) : lset(a_lset), et(a_et), lh(a_lh) {
          D = Dim(et); svs_ptr = make_shared<PointCnt>();
          if(D == 3){//Why is this required?!!
              vector<double> lset_s(lset.Size()); for(int i=0; i<lset.Size(); i++) lset_s[i] = lset[i];
              for(int i=0; i<lset.Size(); i++) lset[i] = lset_s[lset.Size()-1-i];
          }
      }

      void GetIntegrationRule(int order, DOMAIN_TYPE dt, IntegrationRule &intrule);
  };

  const IntegrationRule * StraightCutIntegrationRule(shared_ptr<CoefficientFunction> cf_lset,
                                                     const FlatVector<> & cf_lset_at_element,
                                                     const ElementTransformation & trafo,
                                                     DOMAIN_TYPE dt,
                                                     int intorder,
                                                     LocalHeap & lh);

  class MultiLinearFunction {
  public:
      int D;
      vector<double> c;

      MultiLinearFunction(int a_D=0) : D(a_D) {
          c.resize(1<<D);
          for(int i=0; i<c.size(); i++) c[i] = 0;
      }
      static bool get_bool_i(int h, int a_D, int i) { return (h&(1<<(a_D - i -1))) > 0;}

      static vector<bool> get_bools(int h, int a_D){
          vector<bool> idx(a_D);
          for(int i=0; i<a_D; i++) idx[i] = get_bool_i(h, a_D, i);
          return idx;
      }

      static int get_int(vector<bool> idx) {
          int h=0;
          for(int i=0; i<idx.size(); i++) h += (1<<(idx.size() - i -1))*idx[i];
          return h;
      }

      double& operator[](vector<bool> idx){ return c[get_int(idx)]; }

      template<int Dv>
      double operator()(Vec<Dv> x);

      vector<double> find_root_1D(double x1, double x2);

      MultiLinearFunction get_del_k(int k);

      template<int Dv>
      Vec<Dv> get_grad(Vec<Dv> x);

      void output();

      template<int Dv>
      Vec<2> get_extremal_values_on_hyperrect(Vec<Dv> xL, Vec<Dv> xU);
  };

  template<int D>
  double eval_integrand(Array<MultiLinearFunction> psi, Array<int> s, int k, double x1, double x2, Vec<D-1> x, function<double(Vec<D>)> f, int order);

  template<int D>
  double eval_surface_integrand(MultiLinearFunction phi, int k, double x1, double x2, Vec<D-1> x, function<double(Vec<D>)> f);

}
