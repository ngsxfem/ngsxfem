#include "straightcutrule.hpp"

namespace ngbla {
bool operator==(const Vec<3> a, const Vec<3> b){
    return L2Norm(a - b)<1e-12;
}
}

namespace xintegration
{
  DOMAIN_TYPE CheckIfStraightCut (FlatVector<> cf_lset_at_element, double epsilon) {
    bool haspos = false;
    bool hasneg = false;

    for (auto v : cf_lset_at_element) {
        if (!haspos && (v > epsilon)) haspos = true;
        if (!hasneg && (v < -epsilon)) hasneg = true;
        if(haspos && hasneg) break;
    }

    if ((hasneg && haspos)||(!hasneg && !haspos)) return IF;
    else if (hasneg) return NEG;
    else return POS;
  }

  DOMAIN_TYPE CheckIfStraightCut (vector<double> cf_lset_at_element, double epsilon) {
    bool haspos = false;
    bool hasneg = false;

    for (auto v : cf_lset_at_element) {
        if (!haspos && (v > epsilon)) haspos = true;
        if (!hasneg && (v < -epsilon)) hasneg = true;
        if(haspos && hasneg) break;
    }

    if ((hasneg && haspos)||(!hasneg && !haspos)) return IF;
    else if (hasneg) return NEG;
    else return POS;
  }

  PolytopE SimpleX::CalcIFPolytopEUsingLset(vector<double> lset_on_points){
      //static Timer t ("SimpleX::CalcIFPolytopEUsingLset");
      // RegionTimer reg(t);
      // ThreadRegionTimer reg (t, TaskManager::GetThreadId());
      if(CheckIfStraightCut(lset_on_points) != IF) {
          cout << IM(1) << "Lsetvals: ";
          for(auto d: lset_on_points) cout << d << endl;
          throw Exception ("You tried to cut a simplex with a plain geometry lset function");
      }
      if(D == 1) return SimpleX({Vec<3>(points[0] +(lset_on_points[0]/(lset_on_points[0]-lset_on_points[1]))*(points[1]-points[0]))});
      else {
          Array<Vec<3>> cut_points;
          for(int i = 0; i<points.Size(); i++) {
            for(int j= i+1; j<points.Size(); j++){
                if((lset_on_points[i] >= 0) != (lset_on_points[j] >= 0)){
                    SimpleX s({points[i], points[j]});
                    auto p = s.CalcIFPolytopEUsingLset({lset_on_points[i], lset_on_points[j]});
                    cut_points.Append(p.points[0]);
                }
            }
          }
          return PolytopE(cut_points, D-1);
      }
  }

  double SimpleX::GetVolume(){
      if(D == 0) return 1;
      else if(D == 1) return L2Norm( points[1] - points[0] );
      else if( D == 2) return L2Norm(Cross( Vec<3>(points[2] - points[0]), Vec<3>(points[1] - points[0]) ));
      else if( D == 3) return abs(Determinant<3>(points[3] - points[0], points[2] - points[0], points[1] - points[0]));
      else throw Exception("Calc the Volume of this type of Simplex not implemented!");
  }

  double Quadrilateral::GetVolume(){
      if( D == 2) return L2Norm(Cross( Vec<3>(points[3] - points[0]), Vec<3>(points[1] - points[0])));
      else if( D == 3) return abs(Determinant<3>(points[4] - points[0], points[3] - points[0], points[1] - points[0]));
      else throw Exception("can only handle 2/3 D");
  }

  void SimpleX::GetPlainIntegrationRule(IntegrationRule &intrule, int order) {
      //static Timer t ("SimpleX::GetPlainIntegrationRule");
      // ThreadRegionTimer reg (t, TaskManager::GetThreadId());
      // RegionTimer reg(t);
      double trafofac = GetVolume();

      const IntegrationRule * ir_ngs = nullptr;
      if(D == 0) ir_ngs = & SelectIntegrationRule(ET_POINT, order);
      else if(D == 1) ir_ngs = & SelectIntegrationRule(ET_SEGM, order);
      else if(D == 2) ir_ngs = & SelectIntegrationRule(ET_TRIG, order);
      else if(D == 3) ir_ngs = & SelectIntegrationRule (ET_TET, order);

      for (const auto& ip : *ir_ngs) {
        Vec<3> point(0.0); double originweight = 1.0;
        for (int m = 0; m < points.Size()-1 ;++m) originweight -= ip(m);
          point = originweight * (points[0]);
        for (int m = 0; m < points.Size()-1 ;++m)
          point += ip(m) * (points[m+1]);
        intrule.Append(IntegrationPoint(point, ip.Weight() * trafofac));
      }
  }

  void Quadrilateral::GetPlainIntegrationRule(IntegrationRule &intrule, int order) {
      static Timer t ("Quadrilateral::GetPlainIntegrationRule");
      //RegionTimer reg(t);
      // ThreadRegionTimer reg (t, TaskManager::GetThreadId());
      double trafofac = GetVolume();

      const IntegrationRule * ir_ngs = nullptr;
      Mat<3,3> A(0.0);
      if(D == 2){
          ir_ngs = & SelectIntegrationRule(ET_QUAD, order);
          A(0,0) = points[1][0] - points[0][0]; //TODO: Write down this nicer...
          A(1,0) = points[1][1] - points[0][1];
          A(0,1) = points[3][0] - points[0][0];
          A(1,1) = points[3][1] - points[0][1];
      }
      else if(D == 3){
          ir_ngs =& SelectIntegrationRule(ET_HEX, order);
          //throw Exception("Plain 3D Intrule");
          A.Col(0) = points[1] - points[0];
          A.Col(1) = points[3] - points[0];
          A.Col(2) = points[4] - points[0];
      }

      for(const auto& ip : *ir_ngs){
          Vec<3> point(0.); point = points[0] + A*ip.Point();
          intrule.Append(IntegrationPoint(point, ip.Weight()*trafofac));
      }
  }

  void LevelsetCutSimplex::Decompose(){
      //static Timer t ("LevelsetCutSimplex::Decompose");
      //RegionTimer reg(t);
      // ThreadRegionTimer reg (t, TaskManager::GetThreadId());
      vector<double> lsetvals = lset.initial_coefs;
      PolytopE s_cut = s.CalcIFPolytopEUsingLset(lsetvals);

      if(dt == IF) {
          if(s_cut.points.Size() == s.D) SimplexDecomposition.Append(s_cut);
          else if((s_cut.points.Size() == 4)&&(s.D==3)){
              SimplexDecomposition.Append(SimpleX({s_cut.points[0],s_cut.points[1],s_cut.points[3]}));
              SimplexDecomposition.Append(SimpleX({s_cut.points[0],s_cut.points[2],s_cut.points[3]}));
          }
          else {
            cout << IM(1) << "s.D = " << s.D << " , s_cut.points.Size() = " << s_cut.points.Size() << endl;
            cout << IM(1) << "@ lset vals: " << endl;
            for (auto d: lsetvals) cout << d << endl;
            throw Exception("Bad length of s_cut!");
          }
      }
      else {
          Array<int> relevant_base_simplex_vertices;
          for(int i=0; i<s.D+1; i++)
              if( ((dt == POS) &&(lsetvals[i] >= 0)) || ((dt == NEG) &&(lsetvals[i] < 0)))
                  relevant_base_simplex_vertices.Append(i);
          if((relevant_base_simplex_vertices.Size() == 1)){ //Triangle is cut to a triangle || Tetraeder to a tetraeder
              Array<Vec<3>> point_list(s_cut.points);
              point_list.Append(s.points[relevant_base_simplex_vertices[0]]);
              SimplexDecomposition.Append(SimpleX(point_list));
          }
          else if((relevant_base_simplex_vertices.Size() == 2) && (s.D==2)){ //Triangle is cut to a quad
              Array<Vec<3>> point_listA;
              for(int i : relevant_base_simplex_vertices) point_listA.Append(s.points[i]);
              point_listA.Append(s_cut.points[1]);
              Array<Vec<3>> point_listB(s_cut.points); point_listB.Append(s.points[relevant_base_simplex_vertices[0]]);
              SimplexDecomposition.Append(SimpleX(point_listA));
              SimplexDecomposition.Append(SimpleX(point_listB));
          }
          else if((relevant_base_simplex_vertices.Size() == 2) && (s.D==3)) { //Tetraeder is cut to several tetraeder
              Array<Vec<3>> point_listA(s_cut.points); point_listA[0] = s.points[relevant_base_simplex_vertices[1]];
              Array<Vec<3>> point_listB;
              for(int i : relevant_base_simplex_vertices) point_listB.Append(s.points[i]);
              point_listB.Append(s_cut.points[1]); point_listB.Append(s_cut.points[2]);
              Array<Vec<3>> point_listC(s_cut.points); point_listC[3] = s.points[relevant_base_simplex_vertices[0]];
              SimplexDecomposition.Append(point_listA);
              SimplexDecomposition.Append(point_listB);
              SimplexDecomposition.Append(point_listC);
          }
          else if((relevant_base_simplex_vertices.Size() == 3) && (s.D == 3)){
              Array<Vec<3>> point_listA(s_cut.points); point_listA.Append(s.points[relevant_base_simplex_vertices[2]]);
              Array<Vec<3>> point_listB;
              for(int i : relevant_base_simplex_vertices) point_listB.Append(s.points[i]);
              point_listB.Append(s_cut.points[1]);
              Array<Vec<3>> point_listC(s_cut.points); point_listC[2] = s.points[relevant_base_simplex_vertices[0]];
              point_listC.Append(s.points[relevant_base_simplex_vertices[2]]);
              SimplexDecomposition.Append(point_listA);
              SimplexDecomposition.Append(point_listB);
              SimplexDecomposition.Append(point_listC);
          }
          else {
              cout << IM(1) << "@ lset vals: " << endl;
              for (auto d: lsetvals) cout << d << endl;
              throw Exception("Cutting this part of a tetraeder is not implemented yet!");
          }
      }
  }

  void LevelsetCutSimplex::GetIntegrationRule(IntegrationRule &intrule, int order){
      //static Timer t ("LevelsetCutSimplex::GetIntegrationRule");
      //ThreadRegionTimer reg (t, TaskManager::GetThreadId());
      //RegionTimer reg(t);
      Decompose();
      for(auto s : SimplexDecomposition) s.GetPlainIntegrationRule(intrule, order);
  }

  bool LevelsetCutQuadrilateral::HasTopologyChangeAlongXi(){
      vector<tuple<int,int>> EdgesOfDimXi;
      if( q.D == 2) EdgesOfDimXi = {make_tuple(1,2),make_tuple(0,3)};
      else if (q.D == 3) EdgesOfDimXi = {make_tuple(0,4),make_tuple(1,5),make_tuple(2,6),make_tuple(3,7)};

      Vec<2> vals;
      for (auto t : EdgesOfDimXi) {
          vals[1] = lset(q.points[get<0>(t)]); vals[0] = lset(q.points[get<1>(t)]);
          if(CheckIfStraightCut(vals) == IF) return true;
      }
      return false;
  }

  void LevelsetCutQuadrilateral::Decompose(){
      //static Timer t ("LevelsetCutQuadrilateral::Decompose");
      //RegionTimer reg(t);
      // ThreadRegionTimer reg (t, TaskManager::GetThreadId());
      set<double> TopologyChangeXisS{0,1};
      // int xi = q.D ==2 ? 1 : 2;
      vector<tuple<int,int>> EdgesOfDimXi;
      if( q.D == 2) EdgesOfDimXi = {make_tuple(1,2),make_tuple(0,3)};
      else if (q.D == 3) EdgesOfDimXi = {make_tuple(0,4),make_tuple(1,5),make_tuple(2,6),make_tuple(3,7)};
      else throw Exception("Wrong dimensionality of q in LevelsetCutQuadrilateral::Decompose");
      vector<double> vals(2);
      for (auto t : EdgesOfDimXi) {
          vals[1] = lset(q.points[get<0>(t)]); vals[0] = lset(q.points[get<1>(t)]);
          SimpleX unit_line(ET_SEGM);
          if(CheckIfStraightCut(vals) == IF) {
              TopologyChangeXisS.insert((unit_line.CalcIFPolytopEUsingLset(vals)).points[0][0]);
          }
      }
      vector<double> TopologyChangeXis(TopologyChangeXisS.begin(), TopologyChangeXisS.end());
      sort(TopologyChangeXis.begin(), TopologyChangeXis.end());

      for(int i=0; i<TopologyChangeXis.size() -1; i++){
          double xi0 = TopologyChangeXis[i]; double xi1 = TopologyChangeXis[i+1];
          //if(xi1- xi0 < 1e-12) throw Exception("Orthogonal cut");
          if(q.D == 2){
              array<tuple<double, double>, 2> bnd_vals({make_tuple(q.points[0][0], q.points[2][0]), make_tuple(xi0,xi1)});
              QuadrilateralDecomposition.Append(make_unique<LevelsetCutQuadrilateral>(lset, dt, Quadrilateral(bnd_vals),pol));
          }
          else if(q.D == 3){
              array<tuple<double, double>, 3> bnd_vals({make_tuple(q.points[0][0], q.points[2][0]), make_tuple(q.points[0][1], q.points[2][1]), make_tuple(xi0,xi1)});
              QuadrilateralDecomposition.Append(make_unique<LevelsetCutQuadrilateral>(lset, dt, Quadrilateral(bnd_vals),pol));
          }
      }
  }
  const double c = 0.999;
  const double C = 1./sqrt(1- pow(c,2));

  void LevelsetCutQuadrilateral::GetTensorProductAlongXiIntegrationRule(IntegrationRule &intrule, int order){
      int xi = q.D ==2 ? 1 : 2;

      double xi0 = q.points[0][xi]; double xi1;
      if (xi == 1) xi1 = q.points[2][1];
      else xi1 = q.points[4][2];

      const IntegrationRule & ir_ngs = SelectIntegrationRule(ET_SEGM, order);
      for(auto p1: ir_ngs){
          double xi_ast = xi0 + p1.Point()[0]*(xi1 - xi0);
          IntegrationRule new_intrule;
          vector<double> lsetproj( q.D == 2 ? 2 : 4);
          if(q.D == 2) {
              lsetproj[1] = lset(Vec<3>(q.points[0][0],xi_ast,0)); lsetproj[0] = lset(Vec<3>(q.points[2][0], xi_ast,0));
              LevelsetCutSimplex Codim1ElemAtXast(LevelsetWrapper(lsetproj, ET_SEGM), dt, SimpleX(ET_SEGM));
              Codim1ElemAtXast.GetIntegrationRule(new_intrule, order);
          }
          else if(q.D == 3){
              lsetproj[0] = lset(Vec<3>(q.points[0][0],q.points[0][1], xi_ast)); lsetproj[1] = lset(Vec<3>(q.points[2][0],q.points[0][1], xi_ast));
              lsetproj[2] = lset(Vec<3>(q.points[2][0],q.points[2][1], xi_ast)); lsetproj[3] = lset(Vec<3>(q.points[0][0],q.points[2][1], xi_ast));
              LevelsetCutQuadrilateral Codim1ElemAtXast(LevelsetWrapper(lsetproj, ET_QUAD), dt, Quadrilateral(ET_QUAD), pol);
              Codim1ElemAtXast.GetIntegrationRule(new_intrule, order);
          }
          for(const auto& p2 : new_intrule){
              Vec<3> ip(0.); ip[xi] = xi_ast;
              ip[0] = q.points[0][0] + p2.Point()[0]*(q.points[2][0] - q.points[0][0]);
              if (q.D==3) ip[1] = q.points[0][1] + p2.Point()[1]*(q.points[2][1] - q.points[0][1]);
              double if_scale_factor = 1;
              if (dt == IF){
                  Vec<3> lset_grad(0.); lset_grad = lset.GetGrad(ip);
                  if(q.D == 2) if_scale_factor = L2Norm(lset_grad)/abs(lset_grad[0]);
                  else if(q.D == 3) if_scale_factor = L2Norm(lset_grad)/sqrt(pow(lset_grad[0],2) + pow(lset_grad[1],2));
                  if( isnan(if_scale_factor) || if_scale_factor > C ){
                      cout << IM(1) << "Straightcutrule WARNING: IF scaling factor larger than bound:" << endl;
                      cout << IM(1) << "IF scaling factor: " << if_scale_factor << endl;
                      cout << IM(1) << "dims: " << q.D << endl;
                      cout << IM(1) << "c: " << c << endl;
                      cout << IM(1) << "C: " << C << endl;
                      cout << IM(1) << "This might happen in 3D if the child quad had a bad topology and called its fallback routine." << endl;
                      cout << IM(1) << "If you haven't done so, maybe try POL= OPTIMAL to avoid this" << endl;
                      throw Exception("if_scale_factor larger than bound");
                  }
              }
              intrule.Append(IntegrationPoint( ip , p2.Weight()*p1.Weight()*(xi1-xi0)*if_scale_factor));
          }
      }
  }

  void LevelsetCutQuadrilateral::GetIntegrationRuleOnXYPermutatedQuad(IntegrationRule &intrule, int order){
      IntegrationRule intrule_rotated;
      LevelsetWrapper lset_rotated = lset;
      for(int i : {0,1}) for(int j: {0,1}) for(int k : {0,1}) lset_rotated.c[i][j][k] = lset.c[j][i][k];
      Quadrilateral q_rotated = q;
      for(int i=0; i<q.points.Size(); i++) {
          q_rotated.points[i][0] = q.points[i][1];
          q_rotated.points[i][1] = q.points[i][0];
      }
      Vec<3> tmp = q_rotated.points[1]; q_rotated.points[1] = q_rotated.points[3]; q_rotated.points[3] = tmp;
      if(q.D == 3){ tmp = q_rotated.points[5]; q_rotated.points[5] = q_rotated.points[7]; q_rotated.points[7] = tmp; }
      LevelsetCutQuadrilateral me_rotated(lset_rotated,dt, q_rotated, pol, false);
      me_rotated.GetIntegrationRuleAlongXi(intrule_rotated, order);
      for(const auto& ip: intrule_rotated) intrule.Append(IntegrationPoint(Vec<3>{ip.Point()[1], ip.Point()[0], ip.Point()[2]}, ip.Weight()));
  }


  void LevelsetCutQuadrilateral::GetIntegrationRuleOnXZPermutatedQuad(IntegrationRule &intrule, int order){
      IntegrationRule intrule_rotated;
      LevelsetWrapper lset_rotated = lset;
      for(int i : {0,1}) for(int j: {0,1}) for(int k : {0,1}) lset_rotated.c[i][j][k] = lset.c[k][j][i];
      Quadrilateral q_rotated = q;
      for(int i=0; i<q.points.Size(); i++) {
          q_rotated.points[i][0] = q.points[i][2];
          q_rotated.points[i][2] = q.points[i][0];
      }
      Vec<3> tmp = q_rotated.points[1]; q_rotated.points[1] = q_rotated.points[4]; q_rotated.points[4] = tmp;
      if(q.D == 3) { tmp = q_rotated.points[2]; q_rotated.points[2] = q_rotated.points[7]; q_rotated.points[7] = tmp; }
      LevelsetCutQuadrilateral me_rotated(lset_rotated,dt, q_rotated, pol, false);
      me_rotated.GetIntegrationRule(intrule_rotated, order);
      for(const auto& ip: intrule_rotated) intrule.Append(IntegrationPoint(Vec<3>{ip.Point()[2], ip.Point()[1], ip.Point()[0]}, ip.Weight()));
  }

  void LevelsetCutQuadrilateral::GetIntegrationRuleOnYZPermutatedQuad(IntegrationRule &intrule, int order){
      IntegrationRule intrule_rotated;
      LevelsetWrapper lset_rotated = lset;
      for(int i : {0,1}) for(int j: {0,1}) for(int k : {0,1}) lset_rotated.c[i][j][k] = lset.c[i][k][j];
      Quadrilateral q_rotated = q;
      for(int i=0; i<q.points.Size(); i++) {
          q_rotated.points[i][1] = q.points[i][2];
          q_rotated.points[i][2] = q.points[i][1];
      }
      Vec<3> tmp = q_rotated.points[3]; q_rotated.points[3] = q_rotated.points[4]; q_rotated.points[4] = tmp;
      if(q.D == 3) { tmp = q_rotated.points[2]; q_rotated.points[2] = q_rotated.points[5]; q_rotated.points[5] = tmp; }
      LevelsetCutQuadrilateral me_rotated(lset_rotated,dt, q_rotated, pol, false);
      me_rotated.GetIntegrationRule(intrule_rotated, order);
      for(const auto& ip: intrule_rotated) intrule.Append(IntegrationPoint(Vec<3>{ip.Point()[0], ip.Point()[2], ip.Point()[1]}, ip.Weight()));
  }

  vector<double> LevelsetCutQuadrilateral::GetSufficientCritsQBound(){
      double Vsq = 0;
      vector<Vec<3>> corners = {Vec<3>(0,0,0), Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3>(1,1,0)};
      //auto corners = {Vec<3>(0,0,0), Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3>(1,1,0)};
      vector<int> dim_idx_list = {0,1};
      if(q.D == 3) {
          corners = {Vec<3>(0,0,0), Vec<3>(1,0,0), Vec<3>(0,1,0), Vec<3>(1,1,0),
                     Vec<3>(0,0,1), Vec<3>(1,0,1), Vec<3>(0,1,1), Vec<3>(1,1,1)};
          dim_idx_list = {0,1,2};
      }

      for (auto dim_idx : dim_idx_list) {
          double max = 0;
          for(const auto& p : corners){
              double v = pow(lset.GetGrad(p)[dim_idx], 2);
              if (v > max) max = v;
          }
          Vsq += max;
      }
      double V = sqrt(Vsq);
      vector<double> q_max_of_dim{0,0};
      if(q.D == 3) q_max_of_dim.push_back(0);

      for (const auto& p : corners){
          for (auto dim_idx : dim_idx_list) {
            double q_est = pow(V,2)/(pow(V,2) - pow(lset.GetGrad(p)[dim_idx],2));
            if(q_est > q_max_of_dim[dim_idx]) q_max_of_dim[dim_idx] = q_est;
          }
      }
      for (auto dim_idx : dim_idx_list) q_max_of_dim[dim_idx] =sqrt( 1 - 1/ q_max_of_dim[dim_idx]);
      return q_max_of_dim;
  }

  vector<double> LevelsetCutQuadrilateral::GetExactCritsQBound2D(){
      bool allowance_array[] = {true, true};
      double h_root = -lset.c[1][0][0]/lset.c[1][1][0];
      if ((h_root > 0)&&(h_root < 1)) {
          allowance_array[1] = false;
      }
      h_root = -lset.c[0][1][0]/lset.c[1][1][0];
      if ((h_root > 0)&&(h_root < 1)) {
          allowance_array[0] = false;
      }

      vector<double> q_max_of_dim{0,0};
      for(const auto& p: {Vec<3>{0,0,0}, Vec<3>{1,0,0}, Vec<3>{1,1,0}, Vec<3>{0,1,0}}) {
          auto lset_grad = lset.GetGrad(p);
          double q_y = abs(lset_grad[1])/L2Norm(lset_grad);
          double q_x = abs(lset_grad[0])/L2Norm(lset_grad);
          if (q_y > q_max_of_dim[1]) q_max_of_dim[1] = q_y;
          if (q_x > q_max_of_dim[0]) q_max_of_dim[0] = q_x;
      }
      for (auto idx: {0,1}) if(! allowance_array[idx]) q_max_of_dim[idx] = 2; //Encode non-allowance of a dimension in the double
      return q_max_of_dim;
  }

  DIMENSION_SWAP LevelsetCutQuadrilateral::GetDimensionSwap(){
      if(pol == ALWAYS_NONE) return NONE;
      if(!consider_dim_swap) return ID;

      if(q.D == 2){
        auto Exact_Bound = GetExactCritsQBound2D();

        if(pol == FIRST_ALLOWED){
            if(Exact_Bound[1]< c) return ID;
            else if(Exact_Bound[0] < c) return X_Y;
            else return NONE;
        }
        else if(pol == FIND_OPTIMAL){
            if( (Exact_Bound[0] < c) && (Exact_Bound[1] < c)){
                if(Exact_Bound[1] <= Exact_Bound[0]) return ID;
                else return X_Y;
            }
            else if(Exact_Bound[1] < c) return ID;
            else if(Exact_Bound[0] < c) return X_Y;
            else return NONE;
        }
        else return NONE;        
      }
      else if (q.D == 3){
          auto Suff_Bound = GetSufficientCritsQBound();

          for(auto d : Suff_Bound) if ( isnan(d) ) throw Exception ("Sufficient Criterion calculated nan Bound!");
          if(pol == FIRST_ALLOWED){
              if(Suff_Bound[2] < c) return ID;
              else if(Suff_Bound[1] < c) return Y_Z;
              else if(Suff_Bound[0] < c) return X_Z;
              else return NONE;
          }
          else if(pol == FIND_OPTIMAL){
              int min_dim = distance ( Suff_Bound.begin(), min_element(Suff_Bound.begin(), Suff_Bound.end(), [] (double v1, double v2) {return v1 <= v2;} ));
              if( (min_dim < 0) || (min_dim > 2) ) throw Exception("Finding optimal direction failed");

              if(Suff_Bound[min_dim] < c){
                  if (min_dim == 0) return X_Z;
                  else if(min_dim == 1) return Y_Z;
                  else return ID;
              }
              else return NONE;
          }
          else return NONE;
      }
      else throw Exception("can only handle 2/3 D.");
      return NONE;
  }

  void LevelsetCutQuadrilateral::GetIntegrationRuleAlongXi(IntegrationRule &intrule, int order){
      DOMAIN_TYPE dt_quad = CheckIfStraightCut(q.GetLsetVals(lset));
      if(dt_quad == IF){
          if(HasTopologyChangeAlongXi()) {
              Decompose();
              for(auto& sub_q : QuadrilateralDecomposition){
                  DOMAIN_TYPE dt_decomp_quad = CheckIfStraightCut(sub_q->q.GetLsetVals(lset), 1e-15);
                  if (dt_decomp_quad == IF) sub_q->GetTensorProductAlongXiIntegrationRule(intrule, order);
                  else if (dt_decomp_quad == dt) sub_q->q.GetPlainIntegrationRule(intrule, order);
              }
          }
          else GetTensorProductAlongXiIntegrationRule(intrule, order);
      }
      else if (dt_quad == dt) q.GetPlainIntegrationRule(intrule, order);
  }

  void LevelsetCutQuadrilateral::GetFallbackIntegrationRule(IntegrationRule &intrule, int order){
      vector<vector<int>> sub_simplices;
      if(q.D == 2) sub_simplices = {{0,1,3}, {2,1,3}};
      else if(q.D == 3) sub_simplices = {{3,0,1,5}, {3,1,2,5}, {3,5,2,6}, {4,5,0,3}, {4,7,5,3}, {7,6,5,3}};
      for(auto pnts_idxs: sub_simplices){
          Array<Vec<3>> pnt_list(pnts_idxs.size());
          for(int i=0; i<pnts_idxs.size(); i++) pnt_list[i] = q.points[pnts_idxs[i]];
          SimpleX simpl(pnt_list);
          LevelsetWrapper lset_simpl = lset; lset_simpl.update_initial_coefs(simpl.points);
          DOMAIN_TYPE dt_simpl = CheckIfStraightCut(lset_simpl.initial_coefs);
          if((dt_simpl != IF)&&(dt_simpl == dt)) simpl.GetPlainIntegrationRule(intrule, order);
          else if(dt_simpl == IF) {
              LevelsetCutSimplex trig_cut(lset_simpl, dt, simpl);
              trig_cut.GetIntegrationRule(intrule, order);
          }
      }
  }

  void LevelsetCutQuadrilateral::GetIntegrationRule(IntegrationRule &intrule, int order){
      DIMENSION_SWAP sw = GetDimensionSwap();
      if(sw == ID) GetIntegrationRuleAlongXi(intrule, order);
      else if (sw == X_Y) GetIntegrationRuleOnXYPermutatedQuad(intrule, order);
      else if (sw == Y_Z) GetIntegrationRuleOnYZPermutatedQuad(intrule, order);
      else if (sw == X_Z) GetIntegrationRuleOnXZPermutatedQuad(intrule, order);
      else if (sw == NONE) GetFallbackIntegrationRule(intrule, order);
      else throw Exception ("Unknown Dimension Swap!");
  }

  void LevelsetWrapper::GetCoeffsFromVals(ELEMENT_TYPE et, vector<double> vals){
      Vec<2, Vec<2, Vec<2, double>>> ci;
      for(int i : {0,1}) for(int j: {0,1}) for(int k : {0,1}) ci[i][j][k] = 0.; //TODO: Better Solution??
      if(et == ET_SEGM){
          ci[0][0][0] = vals[1]; ci[1][0][0] = vals[0] - vals[1];
      }
      else if(et == ET_TRIG) {
          ci[0][0][0] = vals[2]; ci[1][0][0] = vals[0] - vals[2]; ci[0][1][0] = vals[1] - vals[2];
      }
      else if(et == ET_TET) {
          ci[0][0][0] = vals[3]; ci[1][0][0] = vals[0] - vals[3]; ci[0][1][0] = vals[1] - vals[3]; ci[0][0][1] = vals[2] - vals[3];
      }
      else if(et == ET_QUAD){
          ci[0][0][0] = vals[0]; ci[1][0][0] = vals[1]-ci[0][0][0], ci[0][1][0] = vals[3] - ci[0][0][0], ci[1][1][0] = vals[2] - ci[1][0][0]- ci[0][1][0]- ci[0][0][0];
      }
      else if(et == ET_HEX){
          ci[0][0][0] = vals[0]; ci[1][0][0] = vals[1]-vals[0]; ci[0][1][0] = vals[3] - vals[0]; ci[0][0][1] = vals[4] - vals[0];
          ci[1][1][0] = vals[2] - ci[1][0][0] - ci[0][1][0] - ci[0][0][0];
          ci[1][0][1] = vals[5] - ci[1][0][0] - ci[0][0][1] - ci[0][0][0];
          ci[0][1][1] = vals[7] - ci[0][1][0] - ci[0][0][1] - ci[0][0][0];
          ci[1][1][1] = vals[6] - ci[1][1][0] - ci[1][0][1] - ci[0][1][1] - ci[1][0][0] - ci[0][1][0] - ci[0][0][1] - ci[0][0][0];
      }
      c = ci; initial_coefs = vals;
  }

  Vec<3> LevelsetWrapper::GetGrad(const Vec<3>& p) const{
      Vec<3> geom_normal(0.);
      geom_normal[0] = c[1][0][0]+c[1][1][0]*p[1]+c[1][0][1]*p[2]+c[1][1][1]*p[1]*p[2];
      geom_normal[1] = c[0][1][0]+c[1][1][0]*p[0]+c[0][1][1]*p[2]+c[1][1][1]*p[0]*p[2];
      geom_normal[2] = c[0][0][1]+c[0][1][1]*p[1]+c[1][0][1]*p[0]+c[1][1][1]*p[0]*p[1];
      return geom_normal;
  }

  Vec<3> LevelsetWrapper::GetNormal(const Vec<3>& p) const{
      Vec<3> geom_normal = GetGrad(p);
      geom_normal /= L2Norm(geom_normal);
      return geom_normal;
  }

  double LevelsetWrapper::operator ()(const Vec<3>& p) const{
      double v = 0;
      for(int i : {0,1}) for(int j: {0,1}) for(int k : {0,1}) v += c[i][j][k]*pow(p[0], i)*pow(p[1], j)*pow(p[2], k);
      return v;
  }

  void LevelsetWrapper::update_initial_coefs(const Array<Vec<3>> &a_points){
      initial_coefs.resize(a_points.Size());
      for(int i=0; i<a_points.Size(); i++){
          double d = operator ()(a_points[i]);
          //initial_coefs[i]= d;
          if(abs(d) > 1e-14) initial_coefs[i]= d;
          else initial_coefs[i] = 1e-14;
      }
  }

  template<unsigned int D>
  void TransformQuadUntrafoToIRInterface(const IntegrationRule & quad_untrafo, const ElementTransformation & trafo, const LevelsetWrapper &lset, IntegrationRule * ir_interface, bool spacetime_mode, double tval){
      for (int i = 0; i < quad_untrafo.Size(); ++i)
      {
          double myweight = quad_untrafo[i].Weight();
          if(spacetime_mode) {
              quad_untrafo[i].SetWeight(tval);
              MarkAsSpaceTimeIntegrationPoint(quad_untrafo[i]);
          }

          MappedIntegrationPoint<D,D> mip(quad_untrafo[i],trafo);
          Mat<D,D> Finv = mip.GetJacobianInverse();

          FlatVector<double> fv(D,&(lset.GetNormal(quad_untrafo[i].Point())(0)));
          Vec<D> normal = Trans(Finv) * fv;
          const double weight = myweight * L2Norm(normal);

          (*ir_interface)[i] = IntegrationPoint (quad_untrafo[i].Point(), weight);
      }
  }

  // integration rules that are returned assume that a scaling with mip.GetMeasure() gives the
  // correct weight on the "physical" domain (note that this is not a natural choice for interface integrals)
  const IntegrationRule * StraightCutIntegrationRule(const FlatVector<> & cf_lset_at_element,
                                                     const ElementTransformation & trafo,
                                                     DOMAIN_TYPE dt,
                                                     int intorder,
                                                     SWAP_DIMENSIONS_POLICY quad_dir_policy,
                                                     LocalHeap & lh,
                                                     bool spacetime_mode,
                                                     double tval)
  {
    static Timer timer("StraightCutIntegrationRule"); 
    RegionTimer reg (timer);

    //static Timer t ("NewStraightCutIntegrationRule");
    // static Timer timercutgeom ("NewStraightCutIntegrationRule::CheckIfCutFast",2);
    // static Timer timermakequadrule("NewStraightCutIntegrationRule::MakeQuadRule",2);

    int DIM = trafo.SpaceDim();

    auto et = trafo.GetElementType();

    if ((et != ET_POINT)&&(et != ET_TRIG)&&(et != ET_TET)&&(et != ET_SEGM)&&(et != ET_QUAD)&&(et != ET_HEX)){
      cout << IM(1) << "Element Type: " << et << endl;
      throw Exception("only trigs, tets, quads for now");
    }
    if ( (et == ET_POINT) && (dt == IF) )
        throw Exception("ET_POINT is only available for volume type ints.");
    bool is_quad = (et == ET_QUAD) || (et == ET_HEX);

    //timercutgeom.Start();
    auto element_domain = CheckIfStraightCut(cf_lset_at_element);
    //timercutgeom.Stop();

    //timermakequadrule.Start();
    IntegrationRule quad_untrafo;
    vector<double> lset_vals(cf_lset_at_element.Size());
    for(int i=0; i<lset_vals.size(); i++) lset_vals[i] = cf_lset_at_element[i];
    LevelsetWrapper lset(lset_vals, et);

    if (element_domain == IF)
    {
      //static Timer timer1("StraightCutElementGeometry::Load+Cut",2);
      //timer1.Start();
      if(!is_quad){
          LevelsetCutSimplex s(lset, dt, SimpleX(et));
          s.GetIntegrationRule(quad_untrafo, intorder);
      }
      else{
          LevelsetCutQuadrilateral q(lset, dt, Quadrilateral(et), quad_dir_policy);
          q.GetIntegrationRule(quad_untrafo, intorder);
          if (quad_untrafo.Size() == 0)
            return nullptr;
      }
      //timer1.Stop();
    }

    const IntegrationRule* ir = nullptr;

    //timermakequadrule.Stop();

    if (element_domain == IF) // there is a cut on the current element
    {
      if (dt == IF)
      {
        auto ir_interface  = new (lh) IntegrationRule(quad_untrafo.Size(),lh);
        if (DIM == 1) TransformQuadUntrafoToIRInterface<1>(quad_untrafo, trafo, lset, ir_interface, spacetime_mode, tval);
        else if (DIM == 2) TransformQuadUntrafoToIRInterface<2>(quad_untrafo, trafo, lset, ir_interface, spacetime_mode, tval);
        else TransformQuadUntrafoToIRInterface<3>(quad_untrafo, trafo, lset, ir_interface, spacetime_mode, tval);
        ir = ir_interface;
      }
      else {
        auto ir_domain = new (lh) IntegrationRule (quad_untrafo.Size(),lh);
        for (int i = 0; i < ir_domain->Size(); ++i)
          (*ir_domain)[i] = IntegrationPoint (quad_untrafo[i].Point(),quad_untrafo[i].Weight());
        ir = ir_domain;
      }
    }
    else
    {
      if (element_domain != dt) //no integration on this element
        return nullptr;
      ir = & (SelectIntegrationRule (trafo.GetElementType(), intorder));
    }

    return ir;
  }

  const IntegrationRule * StraightCutsIntegrationRule(const FlatMatrix<> & cf_lsets_at_element,
                                                     const ElementTransformation & trafo,
                                                     const Array<DOMAIN_TYPE> & dts,
                                                     int intorder,
                                                     SWAP_DIMENSIONS_POLICY quad_dir_policy,
                                                     LocalHeap & lh,
                                                     bool spacetime_mode,
                                                     double tval)
  {

    int DIM = trafo.SpaceDim();

    auto et = trafo.GetElementType();
    int M = cf_lsets_at_element.Width();

    if ((et != ET_TRIG) && (et != ET_TET) && (et != ET_SEGM)){
      cout << IM(1) << "Element Type: " << et << endl;
      throw Exception("only trigs, tets for now");
    }

    IntegrationRule quad_untrafo;

    //This is how to create all the Lset Wrappers ... Collecting them in vector<> or Array<> does not work as there is no empty constructor...
    //Let's see later how we actually will exploit those wrappers...

    auto getLseti_onrefgeom = [&cf_lsets_at_element, &et](int i) {
        vector<double> lset_vals(cf_lsets_at_element.Height());
        for(int ii=0; ii<lset_vals.size(); ii++)
            lset_vals[ii] = cf_lsets_at_element(ii, i);
        LevelsetWrapper lset(lset_vals, et);
        return lset;
    };

    const IntegrationRule* ir = nullptr;

    // outer level
    LevelsetCutSimplex s(getLseti_onrefgeom(0), dts[0], SimpleX(et));
    s.Decompose();
    Array<SimpleX> simplices_at_last_level(s.SimplexDecomposition);

    for (int i = 1; i < M; i++) // all levelset decompositions after the first
    {
      auto lset_on_refgeom = getLseti_onrefgeom(i);
      Array<SimpleX> simplices_at_current_level(0);
      for (auto s : simplices_at_last_level)
      {
        auto lset_on_s = lset_on_refgeom; lset_on_s.update_initial_coefs(s.points);
        auto dt_of_s = CheckIfStraightCut(lset_on_s.initial_coefs);
        if (dt_of_s == IF){
            LevelsetCutSimplex sub_s(lset_on_s, dts[i], s);
            sub_s.Decompose();
            // put sub_s.SimplexDecomposition members to simplices_at_current_level
            simplices_at_current_level.Append(sub_s.SimplexDecomposition);
        }
        else if(dt_of_s == dts[i])
            simplices_at_current_level.Append(s);
      }
      simplices_at_last_level = simplices_at_current_level;
    }

    // Memory leaky:
    // for(auto final_sub_s : simplices_at_last_level)
    //     final_sub_s.GetPlainIntegrationRule(*myir_untrafo, intorder);

    // Workaround:
    IntegrationRule tmp;
    for(auto final_sub_s : simplices_at_last_level)
        final_sub_s.GetPlainIntegrationRule(tmp, intorder);
    
    auto myir_untrafo = new (lh) IntegrationRule(tmp.Size(),lh);
    for (int i = 0; i < tmp.Size(); i++)
      (*myir_untrafo)[i] = tmp[i];
    
    vector<int> dt_is_if_indices;
    for(int i=0; i<M; i++) if ( dts[i] == IF) dt_is_if_indices.push_back(i);

    if(myir_untrafo->Size() == 0) return nullptr;

    if(dt_is_if_indices.size() == 0) //Plain volume case; simple
        ir = myir_untrafo;
    else if(dt_is_if_indices.size() == 1){
        //Do the rescaling according to the one lset function
        auto lset = getLseti_onrefgeom( dt_is_if_indices[0] );

        auto myir = new (lh) IntegrationRule(myir_untrafo->Size(), lh);
        if (DIM == 1) TransformQuadUntrafoToIRInterface<1>(*myir_untrafo, trafo, lset, myir, spacetime_mode, tval);
        else if (DIM == 2) TransformQuadUntrafoToIRInterface<2>(*myir_untrafo, trafo, lset, myir, spacetime_mode, tval);
        else TransformQuadUntrafoToIRInterface<3>(*myir_untrafo, trafo, lset, myir, spacetime_mode, tval);
        ir = myir;
    }
    else if(dt_is_if_indices.size() == 2) {
        // Codim 2 for 2D and 3D
        if(DIM != 2 && DIM != 3) throw Exception("Codim 2 only in 2D and 3D yet!");

        if (DIM == 2){
            for(int i=0; i < myir_untrafo->Size(); i++){
                MappedIntegrationPoint<2,2> mip( (*myir_untrafo)[i],trafo);
                (*myir_untrafo)[i].SetWeight( (*myir_untrafo)[i].Weight() / mip.GetMeasure());
            }
        }
        if(DIM == 3){
            auto lset0 = getLseti_onrefgeom( dt_is_if_indices[0] );
            auto lset1 = getLseti_onrefgeom( dt_is_if_indices[1] );

            auto norm0 = lset0.GetNormal( (*myir_untrafo)[0].Point());
            auto norm1 = lset1.GetNormal( (*myir_untrafo)[0].Point());
            double cp_fac = 1./ L2Norm(Cross(norm0, norm1));

            for(int i=0; i < myir_untrafo->Size(); i++){
                auto old_weight = (*myir_untrafo)[i].Weight();
                MappedIntegrationPoint<3,3> mip( (*myir_untrafo)[i],trafo);

                Mat<3,3> F = mip.GetJacobian();
                Vec<3> normal = cp_fac* F * Cross(norm0, norm1);

                (*myir_untrafo)[i].SetWeight(old_weight*L2Norm(normal) / mip.GetMeasure());
            }
        }
        ir = myir_untrafo;
    }
    else if(dt_is_if_indices.size() == 3){
        // Codim 3 for 3D
        if(DIM != 3) throw Exception("Codim 3 only in 3D!");

        for(int i=0; i < myir_untrafo->Size(); i++){
            MappedIntegrationPoint<3,3> mip( (*myir_untrafo)[i],trafo);
            (*myir_untrafo)[i].SetWeight( (*myir_untrafo)[i].Weight() / mip.GetMeasure());
        }
        ir = myir_untrafo;
    }
    else throw Exception("Codim possibilities available: 2D: 0,1,2; 3D: 0,1,2,3");

    return ir;
  }

  const IntegrationRule * StraightCutIntegrationRuleUntransformed(const FlatVector<> & cf_lset_at_element,
                                                       ELEMENT_TYPE et,
                                                       DOMAIN_TYPE dt,
                                                       int intorder,
                                                       SWAP_DIMENSIONS_POLICY quad_dir_policy,
                                                       LocalHeap & lh)
    {
      static int timer = NgProfiler::CreateTimer ("StraightCutIntegrationRuleUntransformed"); NgProfiler::RegionTimer reg (timer);
      //   static Timer timercutgeom ("NewStraightCutIntegrationRule::CheckIfCutFast");
      //   static Timer timermakequadrule("NewStraightCutIntegrationRule::MakeQuadRule");


      if ((et != ET_POINT)&&(et != ET_TRIG)&&(et != ET_TET)&&(et != ET_SEGM)&&(et != ET_QUAD)&&(et != ET_HEX)){
        cout << IM(1) <<  "Element Type: " << et << endl;
        throw Exception("only trigs, tets, quads for now");
      }
      if ( (et == ET_POINT) && (dt == IF) )
          throw Exception("ET_POINT is only available for volume type ints.");
      bool is_quad = (et == ET_QUAD) || (et == ET_HEX);

      //timercutgeom.Start();
      auto element_domain = CheckIfStraightCut(cf_lset_at_element);
      //timercutgeom.Stop();

      //timermakequadrule.Start();
      IntegrationRule quad_untrafo;
      vector<double> lset_vals(cf_lset_at_element.Size());
      for(int i=0; i<lset_vals.size(); i++) lset_vals[i] = cf_lset_at_element[i];
      LevelsetWrapper lset(lset_vals, et);

      if (element_domain == IF)
      {
        //static Timer timer1("StraightCutElementGeometry::Load+Cut");
        //timer1.Start();
        if(!is_quad){
            LevelsetCutSimplex s(lset, dt, SimpleX(et));
            s.GetIntegrationRule(quad_untrafo, intorder);
        }
        else{
            LevelsetCutQuadrilateral q(lset, dt, Quadrilateral(et), quad_dir_policy);
            q.GetIntegrationRule(quad_untrafo, intorder);
        }
        //timer1.Stop();
      }

      const IntegrationRule* ir = nullptr;

      //timermakequadrule.Stop();

      if (element_domain == IF) // there is a cut on the current element
      {
          auto ir_domain = new (lh) IntegrationRule (quad_untrafo.Size(),lh);
          for (int i = 0; i < ir_domain->Size(); ++i)
            (*ir_domain)[i] = IntegrationPoint (quad_untrafo[i].Point(),quad_untrafo[i].Weight());
          ir = ir_domain;
      }
      else
      {
        if (element_domain != dt) //no integration on this element
          return nullptr;
        ir = & (SelectIntegrationRule (et, intorder));
      }

      return ir;
    }

} // end of namespace
