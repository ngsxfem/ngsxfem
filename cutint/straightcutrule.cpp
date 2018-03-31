#include "straightcutrule.hpp"

namespace ngbla {
bool operator==(const Vec<3> a, const Vec<3> b){
    return L2Norm(a - b)<1e-12;
}
}

namespace xintegration
{
  const bool TRIGGER_MEMORY_LEAK = true;
  const bool SCR_DEBUG_OUTPUT = false; //Temporary solution!!
  const bool SCR_FILE_OUTPUT = false; //Temporary solution!!
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
      static Timer t ("SimpleX::CalcIFPolytopEUsingLset"); RegionTimer reg(t);
      if(CheckIfStraightCut(lset_on_points) != IF) throw Exception ("You tried to cut a simplex with a plain geometry lset function");
      if(D == 1) return SimpleX({Vec<3>(points[0] +(lset_on_points[0]/(lset_on_points[0]-lset_on_points[1]))*(points[1]-points[0]))});
      else {
          vector<Vec<3>> cut_points;
          for(int i = 0; i<points.size(); i++) {
            for(int j= i+1; j<points.size(); j++){
                if((lset_on_points[i] >= 0) != (lset_on_points[j] >= 0)){
                    SimpleX s({points[i], points[j]});
                    auto p = s.CalcIFPolytopEUsingLset({lset_on_points[i], lset_on_points[j]});
                    cut_points.push_back(p.points[0]);
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

  double Quadliteral::GetVolume(){
      if( D == 2) return L2Norm(Cross( Vec<3>(points[3] - points[0]), Vec<3>(points[1] - points[0])));
      else if( D == 3) return abs(Determinant<3>(points[4] - points[0], points[3] - points[0], points[1] - points[0]));
  }

  void SimpleX::GetPlainIntegrationRule(IntegrationRule &intrule, int order) {
      static Timer t ("SimpleX::GetPlainIntegrationRule"); RegionTimer reg(t);
      double trafofac = GetVolume();

      const IntegrationRule * ir_ngs;
      if(D == 0) ir_ngs = & SelectIntegrationRule(ET_POINT, order);
      else if(D == 1) ir_ngs = & SelectIntegrationRule(ET_SEGM, order);
      else if(D == 2) ir_ngs = & SelectIntegrationRule(ET_TRIG, order);
      else if(D == 3) ir_ngs = & SelectIntegrationRule (ET_TET, order);

      for (const auto& ip : *ir_ngs) {
        Vec<3> point(0.0); double originweight = 1.0;
        for (int m = 0; m < points.size()-1 ;++m) originweight -= ip(m);
          point = originweight * (points[0]);
        for (int m = 0; m < points.size()-1 ;++m)
          point += ip(m) * (points[m+1]);
        intrule.Append(IntegrationPoint(point, ip.Weight() * trafofac));
      }
  }

  void Quadliteral::GetPlainIntegrationRule(IntegrationRule &intrule, int order) {
      static Timer t ("Quadliteral::GetPlainIntegrationRule"); RegionTimer reg(t);
      double trafofac = GetVolume();

      const IntegrationRule * ir_ngs;
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

  void LevelsetCuttedSimplex::Decompose(){
      static Timer t ("LevelsetCuttedSimplex::Decompose"); RegionTimer reg(t);
      vector<double> lsetvals = lset.initial_coefs;
      cout << "Simplex decomposition @ lset vals: " << endl;
      for(double d: lsetvals) cout << d << endl;
      PolytopE s_cut = s.CalcIFPolytopEUsingLset(lsetvals);

      if(dt == IF) {
          if(s_cut.points.size() == s.D) SimplexDecomposition.Append(s_cut);
          else if((s_cut.points.size() == 4)&&(s.D==3)){
              SimplexDecomposition.Append(SimpleX({s_cut.points[0],s_cut.points[1],s_cut.points[3]}));
              SimplexDecomposition.Append(SimpleX({s_cut.points[0],s_cut.points[2],s_cut.points[3]}));
          }
          else {
            cout << "s.D = " << s.D << " , s_cut.points.size() = " << s_cut.points.size() << endl;
            cout << "@ lset vals: " << endl;
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
              vector<Vec<3>> point_list(s_cut.points);
              point_list.push_back(s.points[relevant_base_simplex_vertices[0]]);
              SimplexDecomposition.Append(SimpleX(point_list));
          }
          else if((relevant_base_simplex_vertices.Size() == 2) && (s.D==2)){ //Triangle is cut to a quad
              vector<Vec<3>> point_listA;
              for(int i : relevant_base_simplex_vertices) point_listA.push_back(s.points[i]);
              point_listA.push_back(s_cut.points[1]);
              vector<Vec<3>> point_listB(s_cut.points); point_listB.push_back(s.points[relevant_base_simplex_vertices[0]]);
              SimplexDecomposition.Append(SimpleX(point_listA));
              SimplexDecomposition.Append(SimpleX(point_listB));
          }
          else if((relevant_base_simplex_vertices.Size() == 2) && (s.D==3)) { //Tetraeder is cut to several tetraeder
              vector<Vec<3>> point_listA(s_cut.points); point_listA[0] = s.points[relevant_base_simplex_vertices[1]];
              vector<Vec<3>> point_listB;
              for(int i : relevant_base_simplex_vertices) point_listB.push_back(s.points[i]);
              point_listB.push_back(s_cut.points[1]); point_listB.push_back(s_cut.points[2]);
              vector<Vec<3>> point_listC(s_cut.points); point_listC[3] = s.points[relevant_base_simplex_vertices[0]];
              SimplexDecomposition.Append(point_listA);
              SimplexDecomposition.Append(point_listB);
              SimplexDecomposition.Append(point_listC);
          }
          else if((relevant_base_simplex_vertices.Size() == 3) && (s.D == 3)){
              vector<Vec<3>> point_listA(s_cut.points); point_listA.push_back(s.points[relevant_base_simplex_vertices[2]]);
              vector<Vec<3>> point_listB;
              for(int i : relevant_base_simplex_vertices) point_listB.push_back(s.points[i]);
              point_listB.push_back(s_cut.points[1]);
              vector<Vec<3>> point_listC(s_cut.points); point_listC[2] = s.points[relevant_base_simplex_vertices[0]];
              point_listC.push_back(s.points[relevant_base_simplex_vertices[2]]);
              SimplexDecomposition.Append(point_listA);
              SimplexDecomposition.Append(point_listB);
              SimplexDecomposition.Append(point_listC);
          }
          else {
              cout << "@ lset vals: " << endl;
              for (auto d: lsetvals) cout << d << endl;
              throw Exception("Cutting this part of a tetraeder is not implemented yet!");
          }
      }
  }

  void LevelsetCuttedSimplex::GetIntegrationRule(IntegrationRule &intrule, int order){
      static Timer t ("LevelsetCuttedSimplex::GetIntegrationRule"); RegionTimer reg(t);
      Decompose();
      for(auto s : SimplexDecomposition) s.GetPlainIntegrationRule(intrule, order);
  }

  bool LevelsetCuttedQuadliteral::HasTopologyChangeAlongXi(){
      vector<tuple<int,int>> EdgesOfDimXi;
      if( q.D == 2) EdgesOfDimXi = {make_tuple(1,2),make_tuple(0,3)};
      else if (q.D == 3) EdgesOfDimXi = {make_tuple(0,4),make_tuple(1,5),make_tuple(2,6),make_tuple(3,7)};

      Vec<2> vals;
      for (auto t : EdgesOfDimXi) {
          vals[1] = lset(q.points[get<0>(t)]); vals[0] = lset(q.points[get<1>(t)]);
          if(SCR_DEBUG_OUTPUT) {
              cout << "Checking for topology change along edge with lset vals " << vals[0] << " , " << vals[1] << endl;
              cout << "Corresponding points: " << q.points[get<0>(t)] << " \n " << q.points[get<1>(t)] << endl;
          }
          if(CheckIfStraightCut(vals) == IF) return true;
      }
      return false;
  }

  void LevelsetCuttedQuadliteral::Decompose(){
      static Timer t ("LevelsetCuttedQuadliteral::Decompose"); RegionTimer reg(t);
      set<double> TopologyChangeXisS{0,1};
      int xi = q.D ==2 ? 1 : 2;
      vector<tuple<int,int>> EdgesOfDimXi;
      if( q.D == 2) EdgesOfDimXi = {make_tuple(1,2),make_tuple(0,3)};
      else if (q.D == 3) EdgesOfDimXi = {make_tuple(0,4),make_tuple(1,5),make_tuple(2,6),make_tuple(3,7)};
      else throw Exception("Wrong dimensionality of q in LevelsetCuttedQuadliteral::Decompose");
      vector<double> vals(2);
      for (auto t : EdgesOfDimXi) {
          vals[1] = lset(q.points[get<0>(t)]); vals[0] = lset(q.points[get<1>(t)]);
          if(SCR_DEBUG_OUTPUT) cout << "Searching for xi along lset vals " << vals[1] << " , " << vals[0] << endl;
          SimpleX unit_line(ET_SEGM);
          if(CheckIfStraightCut(vals) == IF) {
              TopologyChangeXisS.insert((unit_line.CalcIFPolytopEUsingLset(vals)).points[0][0]);
              if(SCR_DEBUG_OUTPUT) cout << "inserting cut point " << (unit_line.CalcIFPolytopEUsingLset(vals)).points[0][0] << endl;
          }
      }
      vector<double> TopologyChangeXis(TopologyChangeXisS.begin(), TopologyChangeXisS.end());
      sort(TopologyChangeXis.begin(), TopologyChangeXis.end());

      for(int i=0; i<TopologyChangeXis.size() -1; i++){
          double xi0 = TopologyChangeXis[i]; double xi1 = TopologyChangeXis[i+1];
          if(SCR_DEBUG_OUTPUT) cout << "My dim: " << q.D << endl;
          if(SCR_DEBUG_OUTPUT) cout << "Decomposition along interval [xi_0 , xi_1]: " << xi0 << " , " << xi1 << endl;
          //if(xi1- xi0 < 1e-12) throw Exception("Orthogonal cut");
          if(q.D == 2){
              array<tuple<double, double>, 2> bnd_vals({make_tuple(q.points[0][0], q.points[2][0]), make_tuple(xi0,xi1)});
              QuadliteralDecomposition.Append(make_unique<LevelsetCuttedQuadliteral>(lset, dt, Quadliteral(bnd_vals),pol));
          }
          else if(q.D == 3){
              array<tuple<double, double>, 3> bnd_vals({make_tuple(q.points[0][0], q.points[2][0]), make_tuple(q.points[0][1], q.points[2][1]), make_tuple(xi0,xi1)});
              QuadliteralDecomposition.Append(make_unique<LevelsetCuttedQuadliteral>(lset, dt, Quadliteral(bnd_vals),pol));
          }
          if(SCR_DEBUG_OUTPUT) cout << "Child quad successfully created." << endl;
      }
  }
  const double C = 50;
  const double c = sqrt(1-1./pow(C,2));

  void LevelsetCuttedQuadliteral::GetTensorProductAlongXiIntegrationRule(IntegrationRule &intrule, int order){
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
              LevelsetCuttedSimplex Codim1ElemAtXast(LevelsetWrapper(lsetproj, ET_SEGM), dt, SimpleX(ET_SEGM));
              Codim1ElemAtXast.GetIntegrationRule(new_intrule, order);
          }
          else if(q.D == 3){
              lsetproj[0] = lset(Vec<3>(q.points[0][0],q.points[0][1], xi_ast)); lsetproj[1] = lset(Vec<3>(q.points[2][0],q.points[0][1], xi_ast));
              lsetproj[2] = lset(Vec<3>(q.points[2][0],q.points[2][1], xi_ast)); lsetproj[3] = lset(Vec<3>(q.points[0][0],q.points[2][1], xi_ast));
              LevelsetCuttedQuadliteral Codim1ElemAtXast(LevelsetWrapper(lsetproj, ET_QUAD), dt, Quadliteral(ET_QUAD), pol);
              Codim1ElemAtXast.GetIntegrationRule(new_intrule, order);
          }
          for(const auto& p2 : new_intrule){
              Vec<3> ip(0.); ip[xi] = xi_ast;
              ip[0] = q.points[0][0] + p2.Point()[0]*(q.points[2][0] - q.points[0][0]);
              if (q.D==3) ip[1] = q.points[0][1] + p2.Point()[1]*(q.points[2][1] - q.points[0][1]);
              double if_scale_factor = 1;
              if (dt == IF){
                  Vec<3> lset_grad(0.); lset_grad = lset.GetGrad(ip);
                  if(SCR_DEBUG_OUTPUT) cout << "Doing IF scaling: lset_grad: " << lset_grad << endl;
                  if(q.D == 2) if_scale_factor = L2Norm(lset_grad)/abs(lset_grad[0]);
                  else if(q.D == 3) if_scale_factor = L2Norm(lset_grad)/sqrt(pow(lset_grad[0],2) + pow(lset_grad[1],2));
                  if( isnan(if_scale_factor) || if_scale_factor > C ){
                      cout << "IF scaling factor: " << if_scale_factor << endl;
                      throw Exception("if_scale_factor larger than bound");
                  }
                  if(SCR_DEBUG_OUTPUT) cout << "IF scaling factor: " << if_scale_factor << endl;
              }
              intrule.Append(IntegrationPoint( ip , p2.Weight()*p1.Weight()*(xi1-xi0)*if_scale_factor));
          }
      }
  }

  void LevelsetCuttedQuadliteral::GetIntegrationRuleOnXYPermutatedQuad(IntegrationRule &intrule, int order){
      cout << "Permutating X and Y" << endl;
      IntegrationRule intrule_rotated;
      if(SCR_DEBUG_OUTPUT) {
          cout << "Rotation procedure started:" << endl;
          cout << "My lset vals:" << endl;
          for(auto d: q.GetLsetVals(lset)) cout << d << endl;
      }
      LevelsetWrapper lset_rotated = lset;
      for(int i : {0,1}) for(int j: {0,1}) for(int k : {0,1}) lset_rotated.c[i][j][k] = lset.c[j][i][k];
      Quadliteral q_rotated = q;
      for(int i=0; i<q.points.size(); i++) {
          q_rotated.points[i][0] = q.points[i][1];
          q_rotated.points[i][1] = q.points[i][0];
      }
      Vec<3> tmp = q_rotated.points[1]; q_rotated.points[1] = q_rotated.points[3]; q_rotated.points[3] = tmp;
      if(q.D == 3){ tmp = q_rotated.points[5]; q_rotated.points[5] = q_rotated.points[7]; q_rotated.points[7] = tmp; }
      if(SCR_DEBUG_OUTPUT) {
          cout << "rotated lset vals:" << endl;
          for(auto d: q_rotated.GetLsetVals(lset_rotated)) cout << d << endl;
      }
      LevelsetCuttedQuadliteral me_rotated(lset_rotated,dt, q_rotated, pol, false);
      me_rotated.GetIntegrationRuleAlongXi(intrule_rotated, order);
      for(const auto& ip: intrule_rotated) intrule.Append(IntegrationPoint(Vec<3>{ip.Point()[1], ip.Point()[0], ip.Point()[2]}, ip.Weight()));
  }


  void LevelsetCuttedQuadliteral::GetIntegrationRuleOnXZPermutatedQuad(IntegrationRule &intrule, int order){
      cout << "Permutating X and Z" << endl;
      IntegrationRule intrule_rotated;
      if(SCR_DEBUG_OUTPUT) {
          cout << "Rotation procedure started:" << endl;
          cout << "My lset vals:" << endl;
          for(auto d: q.GetLsetVals(lset)) cout << d << endl;
      }
      LevelsetWrapper lset_rotated = lset;
      for(int i : {0,1}) for(int j: {0,1}) for(int k : {0,1}) lset_rotated.c[i][j][k] = lset.c[k][j][i];
      Quadliteral q_rotated = q;
      for(int i=0; i<q.points.size(); i++) {
          q_rotated.points[i][0] = q.points[i][2];
          q_rotated.points[i][2] = q.points[i][0];
      }
      Vec<3> tmp = q_rotated.points[1]; q_rotated.points[1] = q_rotated.points[4]; q_rotated.points[4] = tmp;
      if(q.D == 3) { tmp = q_rotated.points[2]; q_rotated.points[2] = q_rotated.points[7]; q_rotated.points[7] = tmp; }
      if(SCR_DEBUG_OUTPUT) {
          cout << "rotated lset vals:" << endl;
          for(auto d: q_rotated.GetLsetVals(lset_rotated)) cout << d << endl;
      }
      LevelsetCuttedQuadliteral me_rotated(lset_rotated,dt, q_rotated, pol, false);
      me_rotated.GetIntegrationRule(intrule_rotated, order);
      for(const auto& ip: intrule_rotated) intrule.Append(IntegrationPoint(Vec<3>{ip.Point()[2], ip.Point()[1], ip.Point()[0]}, ip.Weight()));
  }

  void LevelsetCuttedQuadliteral::GetIntegrationRuleOnYZPermutatedQuad(IntegrationRule &intrule, int order){
      cout << "Permutating Y and Z" << endl;
      IntegrationRule intrule_rotated;
      if(SCR_DEBUG_OUTPUT) {
          cout << "Rotation procedure started:" << endl;
          cout << "My lset vals:" << endl;
          for(auto d: q.GetLsetVals(lset)) cout << d << endl;
      }
      LevelsetWrapper lset_rotated = lset;
      for(int i : {0,1}) for(int j: {0,1}) for(int k : {0,1}) lset_rotated.c[i][j][k] = lset.c[i][k][j];
      Quadliteral q_rotated = q;
      for(int i=0; i<q.points.size(); i++) {
          q_rotated.points[i][1] = q.points[i][2];
          q_rotated.points[i][2] = q.points[i][1];
      }
      Vec<3> tmp = q_rotated.points[3]; q_rotated.points[3] = q_rotated.points[4]; q_rotated.points[4] = tmp;
      if(q.D == 3) { tmp = q_rotated.points[2]; q_rotated.points[2] = q_rotated.points[5]; q_rotated.points[5] = tmp; }
      if(SCR_DEBUG_OUTPUT) {
          cout << "rotated lset vals:" << endl;
          for(auto d: q_rotated.GetLsetVals(lset_rotated)) cout << d << endl;
      }
      LevelsetCuttedQuadliteral me_rotated(lset_rotated,dt, q_rotated, pol, false);
      me_rotated.GetIntegrationRule(intrule_rotated, order);
      for(const auto& ip: intrule_rotated) intrule.Append(IntegrationPoint(Vec<3>{ip.Point()[0], ip.Point()[2], ip.Point()[1]}, ip.Weight()));
  }

  vector<double> LevelsetCuttedQuadliteral::GetSufficientCritsQBound(){
      double Vsq = 0;
      if(SCR_DEBUG_OUTPUT) cout << "Calculating the suff Crits Q bound in " << q.D << " dims" << endl;
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

  vector<double> LevelsetCuttedQuadliteral::GetExactCritsQBound2D(){
      if(SCR_DEBUG_OUTPUT) cout << "Calculating the exact Crits Q bound in " << q.D << " dims" << endl;
      bool allowance_array[] = {true, true};
      double h_root = -lset.c[1][0][0]/lset.c[1][1][0];
      if ((h_root > 0)&&(h_root < 1)) {
          if(SCR_DEBUG_OUTPUT) cout << "Found the y root " << h_root << " of h!" << endl;
          allowance_array[1] = false;
      }
      h_root = -lset.c[0][1][0]/lset.c[1][1][0];
      if ((h_root > 0)&&(h_root < 1)) {
          if(SCR_DEBUG_OUTPUT) cout << "Found the x root " << h_root << " of h!" << endl;
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

  DIMENSION_SWAP LevelsetCuttedQuadliteral::GetDimensionSwap(){
      if(!TRIGGER_MEMORY_LEAK) {
          cout << "LevelsetCuttedQuadliteral::TransformGeometryIfNecessary on quad\n";
          for(const auto& p : q.points) cout << p << endl;
      }
      //cout << "GetDimensionSwap called in " << q.D << " dimensions with policy " << pol << endl;

      if(pol == ALWAYS_NONE) return NONE;
      if(!consider_dim_swap) return ID;

      if(q.D == 2){
        auto Exact_Bound = GetExactCritsQBound2D();
        //Testing the Sufficient criterion in 2D
        /*cout << "The exact values of q_max are (for both dims) : " << Exact_Bound[0] << "\t" << Exact_Bound[1] << endl;
        auto Suff_Bound = GetSufficientCritsQBound();
        cout << "The sufficient crit. bound is : " << Suff_Bound[0] << "\t" << Suff_Bound[1]  << endl;

        cout << "\t\t !! Ratios between suff and exact: " << Suff_Bound[0] / Exact_Bound[0] << "\t" << Suff_Bound[1] / Exact_Bound[1] << endl;
        */
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
      }
      else if (q.D == 3){
          //return ID;
          auto Suff_Bound = GetSufficientCritsQBound();
          if(SCR_DEBUG_OUTPUT){
            cout << "3D Bounds: " << endl;
            for(auto b: Suff_Bound) cout << b << "\t";
          }
          for(auto d : Suff_Bound) if ( isnan(d) ) throw Exception ("Sufficient Criterion calculated nan Bound!");
          if(pol == FIRST_ALLOWED){
              if(Suff_Bound[2] < c) return ID;
              else if(Suff_Bound[1] < c) return Y_Z;
              else if(Suff_Bound[0] < c) return X_Z;
              else return NONE;
          }
          else if(pol == FIND_OPTIMAL){
              if(SCR_DEBUG_OUTPUT) cout << "Policy: Find_optimal 3D starts ..." << endl;
              int min_dim = distance ( Suff_Bound.begin(), min_element(Suff_Bound.begin(), Suff_Bound.end(), [] (double v1, double v2) {return v1 <= v2;} ));
              if( (min_dim < 0) || (min_dim > 2) ) throw Exception("Finding optimal direction failed");
              if(SCR_DEBUG_OUTPUT) cout << "Min_dim = " << min_dim << " with val : " << Suff_Bound[min_dim] << endl;

              if(Suff_Bound[min_dim] < c){
                  if (min_dim == 0) return X_Z;
                  else if(min_dim == 1) return Y_Z;
                  else return ID;
              }
              else return NONE;
          }
      }
  }

  void LevelsetCuttedQuadliteral::GetIntegrationRuleAlongXi(IntegrationRule &intrule, int order){
      DOMAIN_TYPE dt_quad = CheckIfStraightCut(q.GetLsetVals(lset));
      if(dt_quad == IF){
          if(HasTopologyChangeAlongXi()) {
              Decompose();
              for(auto& sub_q : QuadliteralDecomposition){
                  DOMAIN_TYPE dt_decomp_quad = CheckIfStraightCut(sub_q->q.GetLsetVals(lset), 1e-15);
                  if(SCR_DEBUG_OUTPUT){
                      cout << "Found subquad of type " << dt_decomp_quad << endl;
                      cout << ">>lset vals: " << endl;
                      for(double d: sub_q->q.GetLsetVals(lset)) cout << d << endl;
                  }
                  if (dt_decomp_quad == IF) sub_q->GetTensorProductAlongXiIntegrationRule(intrule, order);
                  else if (dt_decomp_quad == dt) sub_q->q.GetPlainIntegrationRule(intrule, order);
                  cout << "Finished subquad treatment" << endl;
              }
          }
          else GetTensorProductAlongXiIntegrationRule(intrule, order);
      }
      else if (dt_quad == dt) q.GetPlainIntegrationRule(intrule, order);
  }

  void LevelsetCuttedQuadliteral::GetFallbackIntegrationRule(IntegrationRule &intrule, int order){
      vector<vector<int>> sub_simplices;
      if(q.D == 2) sub_simplices = {{0,1,3}, {2,1,3}};
      else if(q.D == 3) sub_simplices = {{3,0,1,5}, {3,1,2,5}, {3,5,2,6}, {4,5,0,3}, {4,7,5,3}, {7,6,5,3}};
      for(auto pnts_idxs: sub_simplices){
          vector<Vec<3>> pnt_list(pnts_idxs.size());
          for(int i=0; i<pnts_idxs.size(); i++) pnt_list[i] = q.points[pnts_idxs[i]];
          SimpleX simpl(pnt_list);
          LevelsetWrapper lset_simpl = lset; lset_simpl.update_initial_coefs(simpl.points);
          DOMAIN_TYPE dt_simpl = CheckIfStraightCut(lset_simpl.initial_coefs);
          if(SCR_DEBUG_OUTPUT){
            cout << "Fallback routine sub levelset vals: " << endl;
            for(auto d: lset_simpl.initial_coefs ) cout << d << endl;
          }
          if((dt_simpl != IF)&&(dt_simpl == dt)) simpl.GetPlainIntegrationRule(intrule, order);
          else if(dt_simpl == IF) {
              LevelsetCuttedSimplex trig_cut(lset_simpl, dt, simpl);
              trig_cut.GetIntegrationRule(intrule, order);
          }
      }
  }

  void LevelsetCuttedQuadliteral::GetIntegrationRule(IntegrationRule &intrule, int order){
      if(!TRIGGER_MEMORY_LEAK) {
          cout << "\n -- LevelsetCuttedQuadliteral::GetIntegrationRule called" << endl;
          cout << "with the lset vals: " << endl;
          //for(auto d: q.GetLsetVals(lset)) cout << d << endl;
          //cout << "on Quad: " << q.points << endl;
          cout << " -- \n" << endl;
      }

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

  void LevelsetWrapper::update_initial_coefs(const vector<Vec<3>> &a_points){
      initial_coefs.resize(a_points.size());
      for(int i=0; i<a_points.size(); i++){
          double d = operator ()(a_points[i]);
          //initial_coefs[i]= d;
          if(abs(d) > 1e-16) initial_coefs[i]= d;
          else initial_coefs[i] = 1e-16;
      }
  }

  template<unsigned int D>
  void TransformQuadUntrafoToIRInterface(const IntegrationRule & quad_untrafo, const ElementTransformation & trafo, const LevelsetWrapper &lset, IntegrationRule * ir_interface){
      for (int i = 0; i < quad_untrafo.Size(); ++i)
      {
          MappedIntegrationPoint<D,D> mip(quad_untrafo[i],trafo);
          Mat<D,D> Finv = mip.GetJacobianInverse();

          Vec<D> normal = Trans(Finv) * lset.GetNormal(quad_untrafo[i].Point());
          const double weight = quad_untrafo[i].Weight() * L2Norm(normal);

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
                                                     LocalHeap & lh)
  {
    static Timer t ("NewStraightCutIntegrationRule");
    static Timer timercutgeom ("NewStraightCutIntegrationRule::CheckIfCutFast");
    static Timer timermakequadrule("NewStraightCutIntegrationRule::MakeQuadRule");

    RegionTimer reg(t);

    int DIM = trafo.SpaceDim();

    auto et = trafo.GetElementType();

    if ((et != ET_TRIG)&&(et != ET_TET)&&(et != ET_SEGM)&&(et != ET_QUAD)&&(et != ET_HEX)){
      cout << "Element Type: " << et << endl;
      throw Exception("only trigs, tets, quads for now");
    }

    bool is_quad = (et == ET_QUAD) || (et == ET_HEX);

    timercutgeom.Start();
    auto element_domain = CheckIfStraightCut(cf_lset_at_element);
    timercutgeom.Stop();

    timermakequadrule.Start();
    IntegrationRule quad_untrafo;
    vector<double> lset_vals(cf_lset_at_element.Size());
    for(int i=0; i<lset_vals.size(); i++) lset_vals[i] = cf_lset_at_element[i];
    LevelsetWrapper lset(lset_vals, et);

    if (element_domain == IF)
    {
      static Timer timer1("StraightCutElementGeometry::Load+Cut");
      timer1.Start();
      if(!is_quad){
          LevelsetCuttedSimplex s(lset, dt, SimpleX(et));
          s.GetIntegrationRule(quad_untrafo, intorder);
      }
      else{
          LevelsetCuttedQuadliteral q(lset, dt, Quadliteral(et), quad_dir_policy);
          q.GetIntegrationRule(quad_untrafo, intorder);
      }
      timer1.Stop();
    }

    const IntegrationRule* ir = nullptr;

    timermakequadrule.Stop();

    if (element_domain == IF) // there is a cut on the current element
    {
      if (dt == IF)
      {
        auto ir_interface  = new (lh) IntegrationRule(quad_untrafo.Size(),lh);
        if (DIM == 1) TransformQuadUntrafoToIRInterface<1>(quad_untrafo, trafo, lset, ir_interface);
        else if (DIM == 2) TransformQuadUntrafoToIRInterface<2>(quad_untrafo, trafo, lset, ir_interface);
        else TransformQuadUntrafoToIRInterface<3>(quad_untrafo, trafo, lset, ir_interface);
        ir = ir_interface;
      }
      else {
        auto ir_domain = new (lh) IntegrationRule (quad_untrafo.Size(),lh);
        for (int i = 0; i < ir_domain->Size(); ++i)
          (*ir_domain)[i] = IntegrationPoint (quad_untrafo[i].Point(),quad_untrafo[i].Weight());
        ir = ir_domain;
      }
      //if(SCR_DEBUG_OUTPUT) cout << "The intrule: " << *ir << endl;
      if(SCR_FILE_OUTPUT){
          ofstream cutrule_outfile("cutrule.dat", fstream::app);
          cutrule_outfile << "dt = " << dt << endl;
          for (auto ip : *ir){
            cutrule_outfile << ip.Point()[0] << "\t" << ip.Point()[1] << "\t" << ip.Point()[2] << "\t\t" << ip.Weight() << endl;
          }
          cutrule_outfile << endl;
          cutrule_outfile.close();
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
} // end of namespace
