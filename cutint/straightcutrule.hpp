#pragma once
#include "xintegration.hpp"
#include <algorithm>
#include <memory>
#include <numeric>

using namespace ngfem;


namespace xintegration
{
  enum SWAP_DIMENSIONS_POLICY {FIRST_ALLOWED, FIND_OPTIMAL};
  enum DIMENSION_SWAP {ID, X_Y, NONE};

  DOMAIN_TYPE CheckIfStraightCut(FlatVector<> cf_lset_at_element, double epsilon = 0);
  DOMAIN_TYPE CheckIfStraightCut(vector<double> cf_lset_at_element, double epsilon = 0);

  class LevelsetWrapper {
  public:
      Vec<2, Vec<2, Vec<2, double>>> c;

      LevelsetWrapper(vector<double> a_vals, ELEMENT_TYPE a_et) { GetCoeffsFromVals(a_et, a_vals); }

      Vec<3> GetNormal(const Vec<3>& p) const;
      Vec<3> GetGrad(const Vec<3>& p) const;
      double operator() (const Vec<3> & p) const;
      vector<double> initial_coefs;
      void update_initial_coefs(const Array<Vec<3>>& a_points);
  private:
      void GetCoeffsFromVals(ELEMENT_TYPE et, vector<double> vals);
  };

  class PolytopE { //The PolytopE which is given as the convex hull of the points
  public:
      Array<Vec<3>> points; //the points
      int D; //Dimension

      PolytopE(const Array<Vec<3>>& a_points, int a_D) : points(a_points), D(a_D) {;}
      PolytopE() { D = -1; } //TODO: Remove this constructor

      vector<double> GetLsetVals(LevelsetWrapper lset){
          vector<double> v;
          for(const auto& p: points) v.push_back(lset(p));
          return v;
      }
  };

  class SimpleX : public PolytopE { //A SimpleX is a PolytopE of dim D with D+1 vertices
  public:
      SimpleX(const Array<Vec<3>>& a_points) : PolytopE(a_points, a_points.Size()-1) {;}
      SimpleX() { D = -1; }

      SimpleX(const PolytopE& p) : PolytopE(p.points, p.D) {
          if(p.D+1 != p.points.Size()) throw Exception("PolytopE -> Simplex constructor called with PolytopE which is not a simplex");
      }

      SimpleX(ELEMENT_TYPE et) {
          if(et == ET_SEGM) *this = SimpleX({{1,0,0}, {0,0,0}});
          else if(et == ET_TRIG) *this = SimpleX({{1,0,0}, {0,1,0}, {0,0,0}});
          else if(et == ET_TET) *this = SimpleX({{1,0,0}, {0,1,0}, {0,0,1}, {0,0,0}});
          else throw Exception ("You tried to create an Simplex with wrong ET");
      }

      PolytopE CalcIFPolytopEUsingLset(vector<double> lset_on_points);

      void GetPlainIntegrationRule(IntegrationRule &intrule, int order);
      double GetVolume();
  };

  class Quadliteral : public PolytopE { //A specific PolytopE: A quadliteral
  public:
      Quadliteral(array<tuple<double, double>, 2> bnds) {
          D = 2; points.SetSize(4);
          points[0] = {get<0>(bnds[0]), get<0>(bnds[1]), 0};
          points[1] = {get<1>(bnds[0]), get<0>(bnds[1]), 0};
          points[2] = {get<1>(bnds[0]), get<1>(bnds[1]), 0};
          points[3] = {get<0>(bnds[0]), get<1>(bnds[1]), 0};
      }
      Quadliteral(array<tuple<double, double>, 3> bnds) {
          D = 3; points.SetSize(8);
          points[0] = {get<0>(bnds[0]), get<0>(bnds[1]), get<0>(bnds[2])};
          points[1] = {get<1>(bnds[0]), get<0>(bnds[1]), get<0>(bnds[2])};
          points[2] = {get<1>(bnds[0]), get<1>(bnds[1]), get<0>(bnds[2])};
          points[3] = {get<0>(bnds[0]), get<1>(bnds[1]), get<0>(bnds[2])};
          points[4] = {get<0>(bnds[0]), get<0>(bnds[1]), get<1>(bnds[2])};
          points[5] = {get<1>(bnds[0]), get<0>(bnds[1]), get<1>(bnds[2])};
          points[6] = {get<1>(bnds[0]), get<1>(bnds[1]), get<1>(bnds[2])};
          points[7] = {get<0>(bnds[0]), get<1>(bnds[1]), get<1>(bnds[2])};
      }
      Quadliteral(ELEMENT_TYPE et) {
          if(et == ET_QUAD) *this = Quadliteral(array<tuple<double, double>, 2>({make_tuple(0,1),make_tuple(0,1)}));
          else if (et == ET_HEX) *this = Quadliteral(array<tuple<double, double>, 3>({make_tuple(0,1),make_tuple(0,1),make_tuple(0,1)}));
          else throw Exception ("You tried to create an Quadliteral with wrong ET");
      }

      void GetPlainIntegrationRule(IntegrationRule &intrule, int order);
      double GetVolume();
  };

  class LevelsetCuttedPolytopE {
  public:
      virtual void GetIntegrationRule(IntegrationRule &intrule, int order) = 0;
      LevelsetWrapper lset;
      DOMAIN_TYPE dt;

      LevelsetCuttedPolytopE(LevelsetWrapper a_lset, DOMAIN_TYPE a_dt): lset(a_lset), dt(a_dt) {;}
  };

  class LevelsetCuttedSimplex : public LevelsetCuttedPolytopE {
  public:
      virtual void GetIntegrationRule(IntegrationRule &intrule, int order);
      SimpleX s;

      LevelsetCuttedSimplex(LevelsetWrapper a_lset, DOMAIN_TYPE a_dt, SimpleX a_s) : LevelsetCuttedPolytopE(a_lset, a_dt), s(a_s) { ;}
  private:
      void Decompose();
      Array<SimpleX> SimplexDecomposition; //TODO: Cf. line 117
  };

  class LevelsetCuttedQuadliteral : public LevelsetCuttedPolytopE {
  public:
      virtual void GetIntegrationRule(IntegrationRule &intrule, int order);
      void GetTensorProductAlongXiIntegrationRule(IntegrationRule &intrule, int order);
      DIMENSION_SWAP GetDimensionSwap(SWAP_DIMENSIONS_POLICY pol = FIRST_ALLOWED);

      Quadliteral q;

      LevelsetCuttedQuadliteral(LevelsetWrapper a_lset, DOMAIN_TYPE a_dt, Quadliteral a_q) : LevelsetCuttedPolytopE(a_lset, a_dt), q(a_q) { ;}
      void GetIntegrationRuleAlongXi(IntegrationRule &intrule, int order);
  private:
      void GetIntegrationRuleOnXYPermutatedQuad(IntegrationRule &intrule, int order);
      void GetIntegrationRuleOnXZPermutatedQuad(IntegrationRule &intrule, int order);
      void GetFallbackIntegrationRule(IntegrationRule &intrule, int order);

      bool HasTopologyChangeAlongXi();
      void Decompose();
      Array<unique_ptr<LevelsetCuttedQuadliteral>> QuadliteralDecomposition;
  };

  template<unsigned int D>
  void TransformQuadUntrafoToIRInterface(const IntegrationRule & quad_untrafo, const ElementTransformation & trafo, const LevelsetWrapper& lset, IntegrationRule * ir_interface);

  const IntegrationRule * StraightCutIntegrationRule(const FlatVector<> & cf_lset_at_element,
                                                     const ElementTransformation & trafo,
                                                     DOMAIN_TYPE dt,
                                                     int intorder,
                                                     LocalHeap & lh);
}
