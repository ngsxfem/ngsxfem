#include "straightcutrule.hpp"

namespace ngbla {
bool operator==(const Vec<3> a, const Vec<3> b){
    return L2Norm(a - b)<1e-12;
}
}

namespace xintegration
{
  DOMAIN_TYPE CheckIfStraightCut (FlatVector<> cf_lset_at_element) {
    bool haspos = false;
    bool hasneg = false;

    for (auto v : cf_lset_at_element) {
        if (!haspos && (v > 1e-10)) haspos = true;
        if (!hasneg && (v < -1e-10)) hasneg = true;
        if(haspos && hasneg) break;
    }

    if (hasneg && haspos) return IF;
    else if (hasneg) return NEG;
    else return POS;
  }

  DOMAIN_TYPE CheckIfStraightCut(const Polytope &s){
      bool haspos = false;
      bool hasneg = false;

      for (int j=0; j<s.Size(); j++) {
          double v = s.GetLset(j);
          if (!haspos && (v > 1e-10)) haspos = true;
          if (!hasneg && (v < -1e-10)) hasneg = true;
          if(haspos && hasneg) break;
      }

      if (hasneg && haspos) return IF;
      else if (hasneg) return NEG;
      else return POS;
  }

  Polytope CalcCutPointLineUsingLset(const Polytope &s){
      if((s.D != 1) ||(s.Size() != 2)) throw Exception("You called the cut-a-line function with a Polytope which is not a line!");

      Vec<3> p = s.GetPoint(0) +(s.GetLset(0)/(s.GetLset(0)-s.GetLset(1)))*(s.GetPoint(1)-s.GetPoint(0));
      s.svs_ptr->Append(make_tuple(p,0));
      return Polytope({s.svs_ptr->Size()-1}, 0, s.svs_ptr);
  }

  Polytope CalcCutPolytopeUsingLset(const Polytope &s){
      Array<int> cut_points;
      if(((s.Size() == 3)&&(s.D == 2))||((s.Size() == 4)&&(s.D == 3))){
        for(int ii = 0; ii<s.Size(); ii++) { int i = s[ii];
            for(int ij= ii+1; ij<s.Size(); ij++){ int j = s[ij];
                if(s.GetLset(ii)*s.GetLset(ij) < -1e-10){
                    Polytope p = CalcCutPointLineUsingLset(Polytope({i, j}, 1, s.svs_ptr));
                    cut_points.Append(p[0]);
                }
            }
        }
    }
    else {
        throw Exception("You tried to cut a Polytope which is not a simplex.");
    }
    return Polytope(cut_points, s.D-1, s.svs_ptr);
  }

  double MeasureSimplVol(const Polytope &s){
      if(s.Size()==2) return L2Norm(Vec<3>(s.GetPoint(1)-s.GetPoint(0)));
      else if(s.Size()==3) return L2Norm(Cross(Vec<3>(s.GetPoint(2) - s.GetPoint(0)), Vec<3>(s.GetPoint(1) - s.GetPoint(0))));
      else if(s.Size()==4) return abs(Determinant<3>(Vec<3>(s.GetPoint(3)-s.GetPoint(0)), Vec<3>(s.GetPoint(2)-s.GetPoint(0)), Vec<3>(s.GetPoint(1)-s.GetPoint(0))));
      else throw Exception("Calc the Volume of this type of Simplex not implemented!");
  }

  void CutSimplexElementGeometry::LoadBaseSimplexFromElementTopology() {
      const POINT3D * verts = ElementTopology::GetVertices(et);

      for(int i=0; i<ElementTopology::GetNVertices(et); i++)
          svs_ptr->Append(make_tuple(Vec<3>{verts[i][0], verts[i][1], verts[i][2]}, lset[i]));

      if((et == ET_TRIG) || (et == ET_TET)){
          Array<int> BaseSimplex(D+1);
          for(int i=0; i<BaseSimplex.Size(); i++) BaseSimplex[i] = i;
          simplices.Append(Polytope(BaseSimplex, D, svs_ptr));
      }
      else throw Exception("Error in LoadBaseSimplexFromElementTopology() - ET_TYPE not supported yet!");
  }

  void CutSimplexElementGeometry::CalcNormal(const Polytope &base_simplex){
      Vec<3> delta_vec;
      double delta_f;
      Vec<3> grad_f; grad_f = 0;
      for(int i=0; i<D; i++) {
          delta_vec = base_simplex.GetPoint(i)-base_simplex.GetPoint(D);
          delta_f = base_simplex.GetLset(i)-base_simplex.GetLset(D);

          for(int j=0; j<3; j++) {
              if (abs(delta_vec[j]) > 1e-10){
                  if(abs(L2Norm(delta_vec) - abs(delta_vec[j]))> 1e-10) throw Exception("Situation to complicated for this type of Grad calculation!");
                  grad_f[j] = delta_f/(delta_vec[j]);
              }
          }
      }
      grad_f /= L2Norm(grad_f);
      normal = grad_f;
  }

  void CutSimplexElementGeometry::CutBaseSimplex(DOMAIN_TYPE dt){
      Polytope s_cut = CalcCutPolytopeUsingLset(simplices[0]);

      if(dt == IF) CalcNormal(simplices[0]);
      simplices.DeleteAll();
      if(dt == IF) {
          if(s_cut.Size() == D) simplices.Append(s_cut);
          else if((s_cut.Size() == 4)&&(D==3)){
              simplices.Append(Polytope({s_cut[0],s_cut[1],s_cut[3]},2, svs_ptr));
              simplices.Append(Polytope({s_cut[0],s_cut[2],s_cut[3]},2, svs_ptr));
          }
          else {
              cout << "s_cut: " << s_cut.ia << endl;
              throw Exception("Bad length of s_cut!");
          }
      }
      else {
          Array<int> nothing;
          Polytope relevant_base_simplex_vertices(nothing, D, svs_ptr);
          for(int i=0; i<D+1; i++)
              if( ((dt == POS) &&(lset[i] > 1e-10)) || ((dt == NEG) &&(lset[i] < -1e-10)))
                  relevant_base_simplex_vertices.Append(i);
          if((relevant_base_simplex_vertices.Size() == 1)){ //Triangle is cut to a triangle || Tetraeder to a tetraeder
              Polytope s(s_cut);
              s.Append(relevant_base_simplex_vertices[0]);
              simplices.Append(s);
          }
          else if((relevant_base_simplex_vertices.Size() == 2) && (D==2)){ //Triangle is cut to a quad
              Polytope s1(relevant_base_simplex_vertices), s2(s_cut);// s1 = ; s2 = ;
              s1.Append(s_cut[1]); s2.Append(relevant_base_simplex_vertices[0]); //The right indices follow from the cutting order
              simplices.Append(s1);
              simplices.Append(s2);
          }
          else if((relevant_base_simplex_vertices.Size() == 2) && (D==3)) { //Tetraeder is cut to several tetraeder
              Polytope s1(s_cut); s1.ia[0] = relevant_base_simplex_vertices[1];
              Polytope s2(relevant_base_simplex_vertices); s2.Append(s_cut[1]); s2.Append(s_cut[2]);
              Polytope s3(s_cut); s3.ia[3] = relevant_base_simplex_vertices[0];
              simplices.Append(s1);
              simplices.Append(s2);
              simplices.Append(s3);
          }
          else if((relevant_base_simplex_vertices.Size() == 3) && (D == 3)){
              Polytope s1(s_cut); s1.Append(relevant_base_simplex_vertices[2]);
              Polytope s2(relevant_base_simplex_vertices); s2.Append(s_cut[1]);
              Polytope s3(s_cut); s3.ia[2] = relevant_base_simplex_vertices[0]; s3.Append(relevant_base_simplex_vertices[2]);
              simplices.Append(s1);
              simplices.Append(s2);
              simplices.Append(s3);
          }
          else {
              throw Exception("Cutting this part of a tetraeder is not implemented yet!");
          }
      }
  }

  void CutSimplexElementGeometry::GetIntegrationRule(int order, DOMAIN_TYPE dt, IntegrationRule &intrule){
      LoadBaseSimplexFromElementTopology();
      CutBaseSimplex(dt);

      int ref_ir_ngs_size; int j =0;
      for(int i=0; i<simplices.Size(); i++) {
        double trafofac = MeasureSimplVol(simplices[i]);

        const IntegrationRule * ir_ngs;
        if(simplices[i].Size() == 2) ir_ngs = & SelectIntegrationRule(ET_SEGM, order);
        else if(simplices[i].Size() == 3) ir_ngs = & SelectIntegrationRule(ET_TRIG, order);
        else if(simplices[i].Size() == 4) ir_ngs = & SelectIntegrationRule (ET_TET, order);

        //cout << "Size of Simplices nr. " << i << ": " << simplices[i].Size() << endl;

        if(i == 0){
            intrule.SetSize(simplices.Size()*ir_ngs->Size());
            ref_ir_ngs_size = ir_ngs->Size();
        }
        else if (ir_ngs->Size() != ref_ir_ngs_size) throw Exception("Different sizes for ir_ngs are not supported!");

        for (auto ip : *ir_ngs) {
          Vec<3> point(0.0); double originweight = 1.0;
          for (int m = 0; m < simplices[i].Size()-1 ;++m) originweight -= ip(m);
            point = originweight * (simplices[i].GetPoint(0));
          for (int m = 0; m < simplices[i].Size()-1 ;++m)
            point += ip(m) * (simplices[i].GetPoint(m+1));
          intrule[j] = IntegrationPoint(point, ip.Weight() * trafofac);
          j++;
        }
      }
  }

  void CutQuadElementGeometry::LoadBaseQuadFromElementTopology() {
      const POINT3D * verts = ElementTopology::GetVertices(et);

      for(int i=0; i<ElementTopology::GetNVertices(et); i++){
          svs_ptr->Append(make_tuple(Vec<3>{verts[i][0], verts[i][1], verts[i][2]}, lset[i]));
          //cout << "Point Nr. " << i << endl;
          //cout << verts[i][0] << "\t" << verts[i][1] << "\t" << verts[i][2] << endl;
      }

      if((et != ET_QUAD)&&(et != ET_HEX)) throw Exception("Error in LoadBaseQuadFromElementTopology() - ET_TYPE not supported yet!");
  }

  Vec<3> CutQuadElementGeometry::GetNormal(const Vec<3>& p) const{
      Vec<3> geom_normal(0.);
      if(D==2){
        geom_normal[0] = lc[1][0][0]+lc[1][1][0]*p[1];
        geom_normal[1] = lc[0][1][0]+lc[1][1][0]*p[0];
      }
      else{
          geom_normal[0] = lc[1][0][0]+lc[1][1][0]*p[1]+lc[1][0][1]*p[2]+lc[1][1][1]*p[1]*p[2];
          geom_normal[1] = lc[0][1][0]+lc[1][1][0]*p[0]+lc[0][1][1]*p[2]+lc[1][1][1]*p[0]*p[2];
          geom_normal[2] = lc[0][0][1]+lc[0][1][1]*p[1]+lc[1][0][1]*p[0]+lc[1][1][1]*p[0]*p[1];
      }
      geom_normal /= L2Norm(geom_normal);
      return geom_normal;
  }

  void CutQuadElementGeometry::FindVolumeAndCutQuads(DOMAIN_TYPE dt){
      function<double(Vec<3>)> levelset = [this] (const Vec<3>& p) {return lc[1][0][0]*p[0]+lc[0][1][0]*p[1]+lc[1][1][0]*p[0]*p[1]+lc[0][0][0];};

      Vec<2, vector<double>> cut_points;
      for (int dim:{0,1}) { cut_points[dim].push_back(0); cut_points[dim].push_back(1);}

      vector<vector<tuple<int,int>>> EdgesOfDim{{make_tuple(0,1),make_tuple(3,2)}, {make_tuple(1,2),make_tuple(0,3)}};
      for(int dim:{0,1}){
          for(tuple<int, int> edge:EdgesOfDim[dim]){
              if(lset[get<0>(edge)]*lset[get<1>(edge)] < -1e-12) {
                  Polytope p = CalcCutPointLineUsingLset(Polytope({get<0>(edge), get<1>(edge)}, 1, svs_ptr));
                  cut_points[dim].push_back(p.GetPoint(0)[dim]);
              }
          }
          sort(cut_points[dim].begin(), cut_points[dim].end());
      }

      for(int i=0; i<cut_points[0].size()-1; i++){
          if(cut_points[0][i+1] -cut_points[0][i] > 1e-10){
              for(int j=0; j<cut_points[1].size()-1; j++){
                  if(cut_points[1][j+1] - cut_points[1][j] > 1e-10){
                      Vec<4, tuple<Vec<3>, double>> quad_points;
                      get<0>(quad_points[0]) = {cut_points[0][i], cut_points[1][j], 0};
                      get<0>(quad_points[1]) = {cut_points[0][i+1], cut_points[1][j], 0};
                      get<0>(quad_points[2]) = {cut_points[0][i+1], cut_points[1][j+1], 0};
                      get<0>(quad_points[3]) = {cut_points[0][i], cut_points[1][j+1], 0};
                      for(int k=0; k<quad_points.Size(); k++) get<1>(quad_points[k]) = levelset(get<0>(quad_points[k]));
                      Polytope poly(quad_points, 2, svs_ptr);
                      DOMAIN_TYPE dt_poly = CheckIfStraightCut(poly);
                      if(dt_poly == IF) Cut_quads.Append(poly);
                      else if(dt_poly == dt) Volume_quads.Append(poly);
                  }
                  else throw Exception("Orthogonal cut y!!");
              }
          }
          else throw Exception("Orthogonal cut x!!");
      }
  }

  void CutQuadElementGeometry::FindVolumeAndCutQuads3D(DOMAIN_TYPE dt){
    function<double(Vec<3>)> levelset = [this] (const Vec<3>& p) {
        double s = 0; for(int i0:{0,1}) for(int i1:{0,1}) for(int i2:{0,1}) s+= lc[i0][i1][i2]*pow(p[0],i0)*pow(p[1],i1)*pow(p[2],i2);
        return s;};

    Vec<3, vector<double>> cut_points;
    for (int dim:{0,1,2}) { cut_points[dim].push_back(0); cut_points[dim].push_back(1);}

    vector<vector<tuple<int,int>>> EdgesOfDim{{make_tuple(0,1),make_tuple(3,2),make_tuple(4,5),make_tuple(7,6)},
                                              {make_tuple(1,2),make_tuple(0,3),make_tuple(4,7),make_tuple(5,6)},
                                              {make_tuple(0,4),make_tuple(1,5),make_tuple(2,6),make_tuple(3,7)}};
    for(int dim:{0,1,2}){
        for(tuple<int, int> edge:EdgesOfDim[dim]){
            if(lset[get<0>(edge)]*lset[get<1>(edge)] < -1e-12) {
                Polytope p = CalcCutPointLineUsingLset(Polytope({get<0>(edge), get<1>(edge)}, 1, svs_ptr));
                cut_points[dim].push_back(p.GetPoint(0)[dim]);
            }
        }
        sort(cut_points[dim].begin(), cut_points[dim].end());
    }
    //cout << "Length of Cut Points x: " << cut_points[0].size() << endl;
    //cout << "Length of Cut Points y: " << cut_points[1].size() << endl;
    //cout << "Length of Cut Points z: " << cut_points[2].size() << endl;

    for(int i=0; i<cut_points[0].size()-1; i++){
        if(cut_points[0][i+1] -cut_points[0][i] > 1e-10){
            for(int j=0; j<cut_points[1].size()-1; j++){
                if(cut_points[1][j+1] - cut_points[1][j] > 1e-10){
                    for(int k=0; k<cut_points[2].size()-1; k++){
                        if(cut_points[2][k+1] - cut_points[2][k] > 1e-10){
                            Vec<8, tuple<Vec<3>, double>> quad_points;
                            get<0>(quad_points[0]) = {cut_points[0][i], cut_points[1][j], cut_points[2][k]};
                            get<0>(quad_points[1]) = {cut_points[0][i+1], cut_points[1][j], cut_points[2][k]};
                            get<0>(quad_points[2]) = {cut_points[0][i+1], cut_points[1][j+1], cut_points[2][k]};
                            get<0>(quad_points[3]) = {cut_points[0][i], cut_points[1][j+1], cut_points[2][k]};
                            get<0>(quad_points[4]) = {cut_points[0][i], cut_points[1][j], cut_points[2][k+1]};
                            get<0>(quad_points[5]) = {cut_points[0][i+1], cut_points[1][j], cut_points[2][k+1]};
                            get<0>(quad_points[6]) = {cut_points[0][i+1], cut_points[1][j+1], cut_points[2][k+1]};
                            get<0>(quad_points[7]) = {cut_points[0][i], cut_points[1][j+1], cut_points[2][k+1]};
                            for(int n=0; n<quad_points.Size(); n++) get<1>(quad_points[n]) = levelset(get<0>(quad_points[n]));
                            Polytope poly(quad_points, 3, svs_ptr);
                            DOMAIN_TYPE dt_poly = CheckIfStraightCut(poly);
                            if(dt_poly == IF) Cut_quads.Append(poly);
                            else if(dt_poly == dt) Volume_quads.Append(poly);
                        }
                        else throw Exception("Orthogonal cut z!!");
                    }
                }
                else throw Exception("Orthogonal cut y!!");
            }
        }
        else throw Exception("Orthogonal cut x!!");
    }
}
  void CutQuadElementGeometry::IntegrateCutQuads(int order, IntegrationRule &intrule) {
      for(auto poly: Cut_quads){
          double x0 = poly.GetPoint(0)[0], x1 = poly.GetPoint(1)[0];
          const IntegrationRule &  ir_ngs = SelectIntegrationRule(ET_SEGM, order);

          for (auto ip : ir_ngs) {
            Vec<3> point(0.0);
            point[0] = x0+ip.Point()[0]*(x1-x0);
            double u = lc[1][0][0]*point[0]+lc[0][0][0], v = lc[1][1][0]*point[0]+lc[0][1][0]; point[1] = -u/v;
            Vec<3> grad_gamma(0.0); grad_gamma[0] = x1-x0;
            grad_gamma[1] = -(x1-x0)*(lc[1][0][0]*v - lc[1][1][0]*u)/(pow(v,2));
            intrule.Append(IntegrationPoint(point, ip.Weight() * L2Norm(grad_gamma)));
          }
      }
  }

  void CutQuadElementGeometry::IntegrateCutQuads3D(int order, IntegrationRule &intrule) {
      for(auto poly: Cut_quads){
          double x0 = poly.GetPoint(0)[0], x1 = poly.GetPoint(2)[0];

          function<double(double)> y_ast = [this](double x) -> double {return -(lc[1][0][0]*x+lc[0][0][0])/(lc[1][1][1]*x+lc[0][1][0]);};
          function<double(double)> Dy_ast = [this](double x) -> double {return -(lc[1][0][0]*(lc[1][1][1]*x+lc[0][1][0]) - lc[1][1][1]*(lc[1][0][0]*x+lc[0][0][0]))/pow(lc[1][1][1]*x+lc[0][1][0], 2);};
          auto y0 = y_ast, y1 = y_ast, Dy0 = Dy_ast, Dy1 = Dy_ast;
          if(poly.GetLset(0) > 1e-12){
              y0 = [&poly] (double x) -> double {return poly.GetPoint(0)[1];};
              Dy0 = [] (double x) -> double {return 0;};
          }
          else {
              y1 = [&poly] (double x) ->double {return poly.GetPoint(2)[1];};
              Dy1 = [] (double x) -> double {return 0;};
          }

          const IntegrationRule & ir_ngs = SelectIntegrationRule(ET_SEGM, order);

          Vec<2> scale_f; scale_f[0] = x1-x0;
          for (auto ip1 : ir_ngs) {
            for(auto ip2: ir_ngs) {
                Vec<3> p(0.0);
                p[0] = x0+ip1.Point()[0]*scale_f[0];
                scale_f[1] = y1(p[0]) - y0(p[0]);
                p[1] = y0(p[0])+ip2.Point()[0]*scale_f[1];
                double u = lc[1][1][0]*p[0]*p[1] +lc[1][0][0]*p[0]+lc[0][1][0]*p[1]+lc[0][0][0];
                double v = lc[1][0][1]*p[0]+lc[0][1][1]*p[1]+lc[0][0][1]+lc[1][0][1]; p[2] = -u/v;
                Vec<3> del_gamma_xi(0.0), del_gamma_eta(0.0);
                del_gamma_xi[0] = scale_f[0];
                del_gamma_eta[1] = scale_f[1];
                del_gamma_xi[1] = ip2.Point()[0]*scale_f[0]*(Dy1(p[0]) - Dy0(p[0]));
                double Dxu = lc[1][0][0]+lc[1][1][0]*p[1]*del_gamma_xi[1]/scale_f[0]+lc[0][1][0]*del_gamma_xi[1]/scale_f[0];
                double Dxv = lc[1][0][1]+lc[1][1][1]*p[1]*del_gamma_xi[1]/scale_f[0]+lc[0][1][1]*del_gamma_xi[1]/scale_f[0];
                del_gamma_xi[2]  = -scale_f[0]*(Dxu*v - Dxv*u)/(pow(v,2));
                del_gamma_eta[2] = -scale_f[1]*((lc[0][1][0] + lc[1][1][0]*p[0])*v - (lc[0][1][1]+lc[1][1][1]*p[0])*u)/(pow(v,2));
                double trafocfac =L2Norm(Cross(del_gamma_xi, del_gamma_eta));
                intrule.Append(IntegrationPoint(p, ip1.Weight() * ip2.Weight()* trafocfac));
                //cout << "Appending the Weight " << ip1.Weight() * ip2.Weight()* trafocfac << " on IF int 3D" << endl;
            }
          }
      }
  }


  void CutQuadElementGeometry::IntegrateVolumeQuads(int order, IntegrationRule &intrule){
      for(auto quad : Volume_quads) {
          Mat<2,2> A(0.0);
          A(0,0) = quad.GetPoint(1)[0] - quad.GetPoint(0)[0];
          A(1,0) = quad.GetPoint(1)[1] - quad.GetPoint(0)[1];
          A(0,1) = quad.GetPoint(3)[0] - quad.GetPoint(0)[0];
          A(1,1) = quad.GetPoint(3)[1] - quad.GetPoint(0)[1];
          const IntegrationRule & ir_ngs = SelectIntegrationRule(ET_QUAD, order);
          double trafofac = abs(Det(A));
          for(auto ip : ir_ngs){
              Vec<3> p(0.); p = quad.GetPoint(0) + A*ip.Point();
              intrule.Append(IntegrationPoint(p, ip.Weight()*trafofac));
          }
      }
  }

  void CutQuadElementGeometry::IntegrateVolumeQuads3D(int order, IntegrationRule &intrule){
      for(auto quad : Volume_quads) {
          Mat<3,3> A(0.0);
          A(0,0) = quad.GetPoint(1)[0] - quad.GetPoint(0)[0];
          A(1,0) = quad.GetPoint(1)[1] - quad.GetPoint(0)[1];
          A(2,0) = quad.GetPoint(1)[2] - quad.GetPoint(0)[2];
          A(0,1) = quad.GetPoint(3)[0] - quad.GetPoint(0)[0];
          A(1,1) = quad.GetPoint(3)[1] - quad.GetPoint(0)[1];
          A(2,1) = quad.GetPoint(3)[2] - quad.GetPoint(0)[2];
          A(0,2) = quad.GetPoint(4)[0] - quad.GetPoint(0)[0];
          A(1,2) = quad.GetPoint(4)[1] - quad.GetPoint(0)[1];
          A(2,2) = quad.GetPoint(4)[2] - quad.GetPoint(0)[2];
          const IntegrationRule & ir_ngs = SelectIntegrationRule(ET_HEX, order);
          double trafofac = abs(Det(A));
          for(auto ip : ir_ngs){
              Vec<3> p(0.); p = quad.GetPoint(0) + A*ip.Point();
              intrule.Append(IntegrationPoint(p, ip.Weight()*trafofac));
              //cout << "Adding the weight " << ip.Weight()*trafofac << " from Volume 3D" << endl;
          }
      }
  }
  /*
  void StraightCutQuadElementGeometry::IntegrateVolumeOfCutQuads3D(DOMAIN_TYPE dt, int order, IntegrationRule &intrule){
      for(auto poly : Cut_quads){
          double x0 = poly.GetPoint(0)[0], x1 = poly.GetPoint(2)[0];
          double y0 = poly.GetPoint(0)[1], y1 = poly.GetPoint(2)[1];

          //function<double(double)> y_ast = [this](double x) -> double {return -(lc[1][0][0]*x+lc[0][0][0])/(lc[1][1][1]*x+lc[0][1][0]);};
          function<double(double, double)> z_ast = [this](double x, double y) -> double {
              return -1.*(lc[1][1][0]*x*y+lc[1][0][0]*x+lc[0][1][0]*y+lc[0][0][0])/(lc[1][0][1]*x+lc[0][1][1]*y+lc[0][0][1]+lc[1][1][1]*x*y);};
          //auto y0 = y_ast, y1 = y_ast;
          cout << "z_ast(0.25,0.25) : " << z_ast(0.25,0.25) << endl;
          function<double(double, double)> z0 = [&poly, &z_ast, &z0] (double x, double y) -> double {return min(poly.GetPoint(4)[2], max(z_ast(x,y), poly.GetPoint(0)[2]));};
          auto z1 = z0;
          if(((dt == POS)&&(poly.GetLset(0) > 1e-12))||((dt == NEG)&&(poly.GetLset(0) < -1e-12))){
              z0 = [&poly] (double x, double y) -> double {return poly.GetPoint(0)[2];};
          }
          else {
              z1 = [&poly] (double x, double y) ->double {return poly.GetPoint(4)[2];};
          }

          IntegrationRule ir_ngs = SelectIntegrationRule(ET_SEGM, order);
          Vec<3> scale_f; scale_f[0] = x1-x0;
          cout << "x interval: " << x0 << " , " << x1 << endl;
          for(auto p1: ir_ngs){
              Vec<3> p(0.); p[0] = x0+p1.Point()[0]*scale_f[0];
              scale_f[1] = y1 - y0;
              cout << "y interval: " << y0 << " , " << y1 << endl;
              for(auto p2 : ir_ngs){
                  p[1] = y0+p2.Point()[0]*scale_f[1];
                  scale_f[2] = z1(p[0], p[1]) - z0(p[0], p[1]);
                  cout << "z interval: " << z0(p[0], p[1]) << " , " << z1(p[0],p[1]) << endl;
                  if(scale_f[2] < 0) throw Exception ("Negative Scale_f!");
                  for(auto p3: ir_ngs){
                      p[2] = z0(p[0],p[1]) + p3.Point()[0]*scale_f[2];
                      intrule.Append(IntegrationPoint(p, p1.Weight()*p2.Weight()*p3.Weight()*scale_f[0]*scale_f[1]*scale_f[2]));
                      cout << "Adding the Weight "<< p1.Weight()*p2.Weight()*p3.Weight()*scale_f[0]*scale_f[1]*scale_f[2] << " from VolumeOfCut3D" << endl;
                  }
              }
          }
      }
  }*/

  void CutQuadElementGeometry::IntegrateVolumeOfCutQuads3D(DOMAIN_TYPE dt, int order, IntegrationRule &intrule){
      for(auto poly : Cut_quads){
          double V_poly = 1; for(int i=0; i<3; i++) V_poly *= poly.GetPoint(6)[i]-poly.GetPoint(0)[i];
          const IntegrationRule & ir_ngs = SelectIntegrationRule(ET_SEGM, order);
          for(auto p1: ir_ngs){
              FlatVector<> lset_proj(4, lh); double xi = p1.Point()[0];
              for(int i=0; i<4; i++) lset_proj[i] = poly.GetLset(i)*xi+poly.GetLset(i+4)*(1-xi);
              //cout << "The levelset Projection: " << lset_proj << endl;
              CutQuadElementGeometry CoDim1Projection(lset_proj, ET_QUAD, lh);
              IntegrationRule ir_tmp;
              CoDim1Projection.GetIntegrationRule(order, dt, ir_tmp);
              for(auto p2: ir_tmp){
                  Vec<3> p(0.);
                  for(int i=0; i<2; i++) p[i] = poly.GetPoint(0)[i]+ p2.Point()[i]*(poly.GetPoint(6)[i]-poly.GetPoint(0)[i]);
                  p[2] = poly.GetPoint(0)[2]+ xi*(poly.GetPoint(6)[2]-poly.GetPoint(0)[2]);
                  intrule.Append(IntegrationPoint(p , p1.Weight()* p2.Weight()* V_poly));
                  //cout << "Appended the weight " << p1.Weight()* p2.Weight()*V_poly << " from VolumeOfCut3D." << endl;
              }
          }
      }
  }

  void CutQuadElementGeometry::IntegrateVolumeOfCutQuads(DOMAIN_TYPE dt, int order, IntegrationRule &intrule){
      for(auto poly : Cut_quads){
          double x0 = poly.GetPoint(0)[0], x1 = poly.GetPoint(2)[0];
          function<double(double)> y_ast = [this](double x) -> double {return -(lc[1][0][0]*x+lc[0][0][0])/(lc[1][1][0]*x+lc[0][1][0]);};
          auto y0 = y_ast, y1 = y_ast;
          if(((dt == POS)&&(poly.GetLset(0)+poly.GetLset(1) > 1e-12))||((dt == NEG)&&(poly.GetLset(0)+poly.GetLset(1) < -1e-12)))
              y0 = [&poly] (double x) -> double {return poly.GetPoint(0)[1];};
          else
              y1 = [&poly] (double x) ->double {return poly.GetPoint(2)[1];};

          const IntegrationRule & ir_ngs = SelectIntegrationRule(ET_SEGM, order);
          Vec<2> scale_f; scale_f[0] = x1-x0;
          for(auto p1: ir_ngs){
              for(auto p2 : ir_ngs){
                  Vec<3> p(0.); p[0] = x0+p1.Point()[0]*scale_f[0];
                  scale_f[1] = y1(p[0]) - y0(p[0]);
                  p[1] = y0(p[0]) + p2.Point()[0]*scale_f[1];
                  intrule.Append(IntegrationPoint(p, p1.Weight()*p2.Weight()*scale_f[0]*scale_f[1]));
                  //cout << "Appended the weight " << p1.Weight()*p2.Weight()*scale_f[0]*scale_f[1] << "from Volume of Cut 2D." << endl;
              }
          }
      }
  }

  void CutQuadElementGeometry::GetIntegrationRule(int order, DOMAIN_TYPE dt, IntegrationRule &intrule){
      LoadBaseQuadFromElementTopology();
      if(D == 2){
          lc[0][0][0] = lset[0]; lc[1][0][0] = lset[1]-lc[0][0][0], lc[0][1][0] = lset[3] - lc[0][0][0], lc[1][1][0] = lset[2] - lc[1][0][0]- lc[0][1][0]- lc[0][0][0];
          for(int i1=0; i1<D; i1++) for(int i2=0; i2<D; i2++) lc[i1][i2][1] = 0.;
      }
      else{
          //cout << "Levelset: " << lset << endl;
          lc[0][0][0] = lset[0]; lc[1][0][0] = lset[1]-lset[0]; lc[0][1][0] = lset[3] - lset[0]; lc[0][0][1] = lset[4] - lset[0];
          lc[1][1][0] = lset[2] - lc[1][0][0] - lc[0][1][0] - lc[0][0][0];
          lc[1][0][1] = lset[5] - lc[1][0][0] - lc[0][0][1] - lc[0][0][0];
          lc[0][1][1] = lset[7] - lc[0][1][0] - lc[0][0][1] - lc[0][0][0];
          lc[1][1][1] = lset[6] - lc[1][1][0] - lc[1][0][1] - lc[0][1][1] - lc[1][0][0] - lc[0][1][0] - lc[0][0][1] - lc[0][0][0];
          //cout << "Lc: " << lc << endl;
      }
      if(D==2) FindVolumeAndCutQuads(dt);
      else FindVolumeAndCutQuads3D(dt);

      if(dt == IF){
          if (D==2) IntegrateCutQuads(order, intrule);
          else IntegrateCutQuads3D(order, intrule);
      }
      else {
          if(D==2) {
              IntegrateVolumeQuads(order, intrule);
              IntegrateVolumeOfCutQuads(dt, order, intrule);
          }
          else{
              IntegrateVolumeQuads3D(order, intrule);
              IntegrateVolumeOfCutQuads3D(dt, order, intrule);
          }
      }
  }

  template<int Dv>
  double MultiLinearFunction::operator()(Vec<Dv> x) {
      if(D != Dv) {
          cout << D << " = D != Dv = " << Dv << endl;
          throw Exception ("x has not the right dim to be parameter of the multilinear function!");
      }
      else {
          double s = 0;
          for(int h=0; h<(1<<D); h++){
              double prod = 1;
              for(int i=0; i<D; i++) prod *= pow(x[i],get_bool_i(h,D,i));
              s += prod*c[h];
          }
          return s;
      }
  }

  vector<double> MultiLinearFunction::find_root_1D(double x1, double x2){
      if(D != 1) return {}; //Exception!
      else if(c[1] == 0) return{};
      else {
          double r = -c[0]/c[1];
          if ((x1<r)&&(r<x2)) return {r};
          else return {};
      }
  }

  MultiLinearFunction MultiLinearFunction::get_del_k(int k) const {
      MultiLinearFunction del_k(D);
      for(int h=0; h<c.size(); h++){
          if(get_bool_i(h,D, k)){
              auto vb = get_bools(h, D);
              vb[k] = false;
              del_k[vb] = c[h];
          }
      }
      return del_k;
  }

  void MultiLinearFunction::output() {
      bool firstoutput = true;
      for(int h=0; h<c.size(); h++){
          if(abs(c[h]) > 1e-12){
              if(!firstoutput) cout << " + ";
              else firstoutput = false;
              cout << c[h];
              for(int i=0; i<D; i++) if (get_bool_i(h,D, i)) cout << "*x"<<i ;
          }
      }
      cout << endl;
  }

  template<int Dv>
  Vec<Dv> MultiLinearFunction::get_grad(Vec<Dv> x) const{
      if(D != Dv) {
          cout << "Get grad: " << D << " = D != Dv = " << Dv << endl;
          throw Exception ("x has not the right dim to be parameter of the multilinear function!");
      }

      Vec<Dv> grad;
      for(int k=0; k<D; k++) {
          grad[k] = get_del_k(k)(x);
      }
      return grad;
  }

  template<int D>
  double eval_integrand(Array<MultiLinearFunction> &psi, Array<int> &s, int k, double x1, double x2, Vec<D-1> x, function<double(Vec<D>)> f, int order) {
      vector<double> R{x1,x2};
      for(auto psi_i : psi){
          MultiLinearFunction psi_i_new(1);
          for(int h=0; h<psi_i.c.size(); h++) {
              vector<bool> idx = MultiLinearFunction::get_bools(h,D);
              double prod = 1;
              for(int j=0; j<k; j++) prod *= pow(x[j],idx[j]);
              for(int j=k; j<D-1; j++) prod *= pow(x[j],idx[j+1]);
              psi_i_new.c[idx[k]] += psi_i[idx]*prod;
          }
          cout << "psi_i_new: " << endl; psi_i_new.output();
          auto rv = psi_i_new.find_root_1D(x1,x2);
          R.insert(R.begin(), rv.begin(), rv.end());
      }
      sort(R.begin(), R.end());
      double I =0;
      for(int j=0; j<R.size()-1; j++){
          double L = R[j+1]-R[j]; Vec<D> xc;
          for(int i=0; i<k; i++) xc[i] = x[i];
          xc[k] = 0.5*(R[j+1]+R[j]);
          for(int i=k+1; i<D; i++) xc[i] = x[i-1];
          bool cond = true;
          for(int i=0; i<psi.Size(); i++) if(s[i]*psi[i](xc) < 0) cond = false;
          if(cond){
              const IntegrationRule & ir_ngs = SelectIntegrationRule(ET_SEGM, order);
              for(auto ip: ir_ngs){
                  Vec<D> p = xc; p[k] = R[j] + L*ip.Point()[0];
                  I += L*ip.Weight()*f(p);
              }
          }
      }
      return I;
  }

  template<int D>
  void eval_integrand(Array<MultiLinearFunction> &psi, Array<int> &s, int k, double x1, double x2, Vec<D-1> x, int order, IntegrationRule& result) {
      vector<double> R{x1,x2};
      for(auto psi_i : psi){
          MultiLinearFunction psi_i_new(1);
          for(int h=0; h<psi_i.c.size(); h++) {
              vector<bool> idx = MultiLinearFunction::get_bools(h,D);
              double prod = 1;
              for(int j=0; j<k; j++) prod *= pow(x[j],idx[j]);
              for(int j=k; j<D-1; j++) prod *= pow(x[j],idx[j+1]);
              psi_i_new.c[idx[k]] += psi_i[idx]*prod;
          }
          cout << "eval_integrand psi_i_new: " << endl; psi_i_new.output();
          auto rv = psi_i_new.find_root_1D(x1,x2);
          R.insert(R.begin(), rv.begin(), rv.end());
      }
      sort(R.begin(), R.end());
      for(int j=0; j<R.size()-1; j++){
          double L = R[j+1]-R[j]; Vec<D> xc;
          for(int i=0; i<k; i++) xc[i] = x[i];
          xc[k] = 0.5*(R[j+1]+R[j]);
          for(int i=k+1; i<D; i++) xc[i] = x[i-1];
          bool cond = true;
          for(int i=0; i<psi.Size(); i++) if(s[i]*psi[i](xc) < 0) cond = false;
          if(cond){
              const IntegrationRule & ir_ngs = SelectIntegrationRule(ET_SEGM, order);
              for(auto ip: ir_ngs){
                  Vec<D> p = xc; p[k] = R[j] + L*ip.Point()[0];
                  IntegrationPoint ip_new(p, L*ip.Weight());
                  result.Append(ip_new);
              }
          }
      }
  }

  template<int D>
  double eval_surface_integrand(MultiLinearFunction phi, int k, double x1, double x2, Vec<D-1> x, function<double(Vec<D>)> f) {
      MultiLinearFunction psi_new(1);
      for(int h=0; h<phi.c.size(); h++) {
          vector<bool> idx = MultiLinearFunction::get_bools(h, D);
          //for(int j=0; j<D; j++) idx[j] = idx[j];
          double prod = 1;
          for(int j=0; j<k; j++) prod *= pow(x[j],idx[j]);
          for(int j=k; j<D-1; j++) prod *= pow(x[j],idx[j+1]);
          psi_new.c[idx[k]] += phi[idx]*prod;
      }
      cout << "Psi new eval_surface_integrand: " << endl;
      psi_new.output();
      auto rv = psi_new.find_root_1D(x1,x2);
      if(rv.size() == 0) return 0;
      else {
          Vec<D> p;
          for(int i=0; i<k; i++) p[i] = x[i];
          p[k] = rv[0];
          for(int i=k+1; i<D; i++) p[i] = x[i-1];
          return f(p)*L2Norm(phi.get_grad(p))/abs(phi.get_del_k(k)(p));
      }
  }

  template<int D>
  void eval_surface_integrand(MultiLinearFunction phi, int k, double x1, double x2, Vec<D-1> x, IntegrationRule& result) {
      MultiLinearFunction psi_new(1);
      for(int h=0; h<phi.c.size(); h++) {
          vector<bool> idx = MultiLinearFunction::get_bools(h, D);
          double prod = 1;
          for(int j=0; j<k; j++) prod *= pow(x[j],idx[j]);
          for(int j=k; j<D-1; j++) prod *= pow(x[j],idx[j+1]);
          psi_new.c[idx[k]] += phi[idx]*prod;
      }
      cout << "eval_surface_integrand psi_new: " << endl;
      psi_new.output();
      auto rv = psi_new.find_root_1D(x1,x2);
      if(rv.size() > 0){
          Vec<D> p;
          for(int i=0; i<k; i++) p[i] = x[i];
          p[k] = rv[0];
          for(int i=k+1; i<D; i++) p[i] = x[i-1];
          result.Append(IntegrationPoint(p, L2Norm(phi.get_grad(p))/abs(phi.get_del_k(k)(p))));
      }
  }

  template<int Dv>
  double MultiLinearFunction::get_largest_abs_on_hyperrect(Vec<Dv> xL, Vec<Dv> xU){
      if(Dv != D) throw Exception ("Dimension mismatch!");
      vector<double> vals(1<<Dv);
      for(int h=0; h<vals.size(); h++) {
          Vec<Dv> p;
          for(int i=0; i<D; i++){
              if(MultiLinearFunction::get_bool_i(h,Dv,i)) p[i] = xU[i];
              else p[i] = xL[i];
          }
          vals[h] = abs(operator ()(p));
      }
      auto res = max_element(vals.begin(), vals.end());
      return *res;
  }

  template<int Dv>
  double PolynomeFunction::operator()(Vec<Dv> x){
      if (Dv != D) throw Exception("Dimension mismatch!");
      double val = 0;
      for(auto c_tuple: c){
          auto exponents = get<0>(c_tuple);
          double coeff = get<1>(c_tuple);

          double prod = 1;
          for(int i=0; i<D; i++) prod *= pow( x[i], exponents[i]);
          val += coeff * prod;
      }
      return val;
  }

  vector<double> PolynomeFunction::find_root_1D(double x1, double x2, int subdivs=50, int bisection_iterations = 40){
      vector<double> vals(subdivs); vector<tuple<double,double>> sign_change_intervals;
      vector<double> roots; double delta_x =  (x2 - x1)/subdivs;
      for(int i=0; i<subdivs; i++){
          double xi = x1 + delta_x*i;
          vals[i] = operator()(Vec<1>{xi});
          if(vals[i] == 0) roots.push_back(xi);
          if(i >= 1) if(vals[i-1] * vals[i] < 0) sign_change_intervals.push_back(make_tuple( xi-delta_x, xi));
      }
      for(auto interval : sign_change_intervals){
          double a = get<0>(interval), b = get<1>(interval); double x_mid;
          double aval = operator()(Vec<1>{a}), bval = operator()(Vec<1>{b});;
          for(int j=0; j<bisection_iterations; j++){
              x_mid = 0.5*(a+b);
              double val = operator() (Vec<1>{x_mid});
              if(val == 0) break;
              if(val * aval < 0) {
                  b = x_mid; bval = val;
              }
              else if(val * bval < 0){
                  a = x_mid; aval = val;
              }
              else throw Exception("Strange sign structure during bisection!");
          }
          roots.push_back(0.5*(a+b));
      }
      return roots;
  }

  PolynomeFunction PolynomeFunction::get_del_k(int k) const{
      PolynomeFunction del_k(D);
      for(auto c_tuple : c){
          auto exponents = c_tuple.first;
          double coeff = c_tuple.second;

          if(exponents[k] >= 1){
              exponents[k] -= 1;
              del_k.c[exponents] = coeff*(exponents[k] +1); //.insert(make_tuple( , ));
          }
      }
      return del_k;
  }

  void PolynomeFunction::output(){
      bool firstoutput = true;
      for(auto c_tuple : c){
          if(!firstoutput) cout << " + ";
          else firstoutput = false;

          auto exponents = c_tuple.first;
          double coeff = c_tuple.second;
          cout << coeff;
          for(int i=0; i<exponents.size(); i++) cout << "*x" << i << "^" << exponents[i];
      }
      cout << endl;
  }

  template<int Dv>
  Vec<Dv> PolynomeFunction::get_grad(Vec<Dv> x) const {
      if(D != Dv) throw Exception ("x has not the right dim to be parameter of the multilinear function!");

      Vec<Dv> grad;
      for(int k=0; k<D; k++) {
          grad[k] = get_del_k(k)(x);
      }
      return grad;
  }

  void DebugPolynomeClass(){
      PolynomeFunction p(2);
      p.c[vector<int>{0,0}] = 0.53;
      p.c[vector<int>{9,1}] = 0.15;

      cout << "c(0.3,0.4): " << p(Vec<2>{0.3,0.4}) << endl;

      PolynomeFunction p2(1);
      p2.c[vector<int>{0}] = -2.;
      p2.c[vector<int>{2}] = +1.;
      cout << "sqrt(2) = ";
      auto sqrt2 = p2.find_root_1D(0,2);
      for (auto d: sqrt2) cout << d << endl;

      PolynomeFunction p3(1);
      p3.c[vector<int>{0}] = -20.;
      p3.c[vector<int>{3}] = -30.;
      p3.c[vector<int>{4}] = 5.;
      cout << "roots of 5*x^4 - 30*x^3 - 20: ";
      auto roots = p3.find_root_1D(-10,10);
      for (auto d: roots) cout << d << endl;

      cout << "p3: "; p3.output();
      cout << "p3 ': "; p3.get_del_k(0).output();
  }

  template<class T>
  inline constexpr T pow(const T base, unsigned const exponent)
  {
      return (exponent == 0) ? 1 : (base * pow(base, exponent-1));
  }

  template<int D>
  void Get_Tensor_Product_IR(int order, Vec<D> xL, Vec<D> xU, IntegrationRule& result){
      const IntegrationRule & ir_ngs = SelectIntegrationRule(ET_SEGM, order);
      int q = ir_ngs.Size();
      for(int h=0; h<pow(q, D); h++){
          //for(int i=0; i<D; i++) cout << (h%pow(q,i+1))/pow(q,i);
          //cout << endl;
          double w = 1; Vec<D> point;
          for(int i=0; i<D; i++){
              auto ip = ir_ngs[(h%pow(q,i+1))/pow(q,i)];
              w *= ip.Weight();
              point[i] = xL[i] + (xU[i]-xL[i])*ip.Point()[0];
          }
          result.Append(IntegrationPoint(point, w));
      }
  }

  int sgn(int m, int s, bool S, int sigma){
      if( (m == sigma*s) || S) return sigma*m;
      else return 0;
  }

  auto sgn_L(int m, int s, bool S) {return sgn(m,s,S,-1);}
  auto sgn_U(int m, int s, bool S) {return sgn(m,s,S,+1);}

  template<int D>
  double integrate_saye(Array<MultiLinearFunction>& psi, Array<int>& s, Vec<D> xL, Vec<D> xU, function<double(Vec<D>)> f, bool S, int order, int subdivlevel = 0) {
    Vec<D> xc; for(int i=0; i<D; i++) xc[i] = 0.5*(xL[i]+xU[i]);
    Array<MultiLinearFunction> psi_pruned; Array<int> s_pruned;
    for(int i=psi.Size()-1; i>=0; i--){
        auto psi_c = psi[i]; psi_c.c[0] -= psi[i](xc);
        auto delta = psi_c.get_largest_abs_on_hyperrect(xL, xU);
        if(abs(psi[i](xc)) >= delta){
            if(s[i]*psi[i](xc) < 0) return 0.;
        }
        else {
            psi_pruned.Append(psi[i]); s_pruned.Append(s[i]);
        }
    }
    if(psi_pruned.Size() == 0){
        double vol_U = 1;
        for(int i=0; i<D; i++) vol_U *= (xU[i] - xL[i]);
        double I = 0;
        IntegrationRule ir;
        Get_Tensor_Product_IR(order, xL, xU, ir);
        for(auto ip:ir) I += vol_U*ip.Weight()*f(ip.Point());
        return I;
    }
    vector<double> partial_derivs(D);
    for(int i=0; i<D; i++) partial_derivs[i] = abs(psi_pruned[0].get_del_k(i)(xc));
    int k = distance(partial_derivs.begin(), max_element(partial_derivs.begin(), partial_derivs.end()));
    Array<MultiLinearFunction> psitilde; Array<int> stilde;
    for(int i=0; i<psi_pruned.Size(); i++){
        auto psi_i = psi_pruned[i];
        Vec<D> g = psi_i.get_grad(xc);
        Vec<D> delta;
        for(int j=0; j<D; j++){
            auto psi_c = psi_i.get_del_k(j); psi_c.c[0] -= g[j];
            delta[j] = psi_c.get_largest_abs_on_hyperrect(xL,xU);
        }
        double sum=0; for(int j=0; j<D; j++) sum += pow(g[j]+delta[j],2);
        if( (abs(g[k]) > delta[k]) && (sum / pow(g[k]-delta[k],2) < 20.)){
            MultiLinearFunction psi_i_L(D-1), psi_i_U(D-1);
            for(int h = 0; h<psi_i_L.c.size(); h++){
                auto bools = MultiLinearFunction::get_bools(h, D-1);
                bools.insert(bools.begin()+k, 0);
                psi_i_L.c[h] += psi_i[bools]; psi_i_U.c[h] += psi_i[bools];
                bools = MultiLinearFunction::get_bools(h, D-1);
                bools.insert(bools.begin()+k, 1);
                psi_i_L.c[h] += psi_i[bools]*xL[k]; psi_i_U.c[h] += psi_i[bools]*xU[k];
            }
            psitilde.Append(psi_i_L); psitilde.Append(psi_i_U);
            int sign_gk = 0;
            if(g[k] < 0) sign_gk = -1;
            else if(g[k] > 0) sign_gk = 1;
            stilde.Append(sgn_L(sign_gk, s_pruned[i], S)); stilde.Append(sgn_U(sign_gk, s_pruned[i], S));
        }
        else {
            if(subdivlevel >= 10) throw Exception ("Saye asked for Subdivlevel 11!");
            vector<double> d(D);
            for(int j=0; j<D; j++) d[j] = abs(xU[j] - xL[j]);
            int k = distance(d.begin(), max_element(d.begin(), d.end()));
            cout << "Dimension for subdivision: " << k << endl;
            Vec<D> xL1, xL2, xU1, xU2;
            for(int j=0; j<D; j++){
                if(j == k){
                    xL1[j] = xL[j]; xU1[j] = 0.5*(xL[j] + xU[j]);
                    xL2[j] = xU1[j]; xU2[j] = xU[j];
                }
                else{
                    xL1[j] = xL[j]; xL2[j] = xL[j];
                    xU1[j] = xU[j]; xU2[j] = xU[j];
                }
            }
            cout << "xL1: " << xL1 << endl << "xU1: " << xU1 << endl << "xL2: " << xL2 << endl << "xU2: " << xU2 << endl;
            return integrate_saye(psi_pruned, s_pruned, xL1, xU1, f, S, order, subdivlevel+1) + integrate_saye(psi_pruned, s_pruned, xL2, xU2, f, S, order, subdivlevel+1);
        }
    }
    function<double(Vec<D-1>)> ftilde;
    if(S) ftilde = [&] (Vec<D-1> x) {return eval_surface_integrand<D>(psi_pruned[0],k, xL[k], xU[k], x, f); };
    else  ftilde = [&] (Vec<D-1> x) {return eval_integrand<D>(psi_pruned, s_pruned, k, xL[k], xU[k], x, f, order); };

    Vec<D-1> xLtilde, xUtilde;
    for(int i=0; i<D; i++){
        if(i < k) { xLtilde[i] = xL[i]; xUtilde[i] = xU[i]; }
        else if (i > k) { xLtilde[i-1] = xL[i]; xUtilde[i-1] = xU[i]; };
    }

    return integrate_saye<D-1>(psitilde, stilde, xLtilde, xUtilde, ftilde, false, order);
  }
  template<>
  double integrate_saye<1>(Array<MultiLinearFunction>& psi, Array<int>& s, Vec<1> xL, Vec<1> xU, function<double(Vec<1>)> f, bool S, int order, int subdivlevel) {
    return eval_integrand<1>(psi, s, 0, xL[0], xU[0], {}, f, order);
  }

  template<int D>
  void integrate_saye(Array<MultiLinearFunction>& psi, Array<int>& s, Vec<D> xL, Vec<D> xU, bool S, int order, IntegrationRule& result, int subdivlevel = 0) {
    if(subdivlevel > 0) cout << "Starting Saye with Subdivlevel " << subdivlevel << endl;
    Vec<D> xc; for(int i=0; i<D; i++) xc[i] = 0.5*(xL[i]+xU[i]);
    Array<MultiLinearFunction> psi_pruned; Array<int> s_pruned;
    for(int i=psi.Size()-1; i>=0; i--){
        auto psi_c = psi[i]; psi_c.c[0] -= psi[i](xc);
        auto delta = psi_c.get_largest_abs_on_hyperrect(xL, xU);
        if(abs(psi[i](xc)) >= delta){
            if(s[i]*psi[i](xc) < 0) return;
        }
        else {
            psi_pruned.Append(psi[i]); s_pruned.Append(s[i]);
        }
    }
    double vol_U = 1;
    for(int i=0; i<D; i++) vol_U *= (xU[i] - xL[i]);
    if(psi_pruned.Size() == 0){
        IntegrationRule ir;
        Get_Tensor_Product_IR(order, xL, xU, ir);
        for(auto ip:ir) result.Append(IntegrationPoint(ip.Point(), ip.Weight()*vol_U));
        return;
    }
    vector<double> partial_derivs(D);
    for(int i=0; i<D; i++) {
        partial_derivs[i] = abs(psi_pruned[0].get_del_k(i)(xc));
        //cout << "Partial Deriv " << i << " : " << partial_derivs[i] << endl;
    }
    int k = distance(partial_derivs.begin(), max_element(partial_derivs.begin(), partial_derivs.end()));
    cout << "k: " << k << endl;
    Array<MultiLinearFunction> psitilde; Array<int> stilde;
    for(int i=0; i<psi_pruned.Size(); i++){
        auto psi_i = psi_pruned[i];
        Vec<D> g = psi_i.get_grad(xc);
        Vec<D> delta;
        for(int j=0; j<D; j++){
            MultiLinearFunction psi_c(D); psi_c = psi_i.get_del_k(j); psi_c.c[0] -= g[j];
            delta[j] = psi_c.get_largest_abs_on_hyperrect(xL,xU);
        }
        cout << "Delta: " << delta << endl; cout << "g: " << g << endl;
        double sum=0; for(int j=0; j<D; j++) sum += pow(g[j]+delta[j],2);
        //cout << "sum: " << sum << endl;
        //cout << " sum / pow(g[k]-delta[k],2) " << sum / pow(g[k]-delta[k],2.) << endl;
        if( (abs(g[k]) > delta[k]) && (sum / pow(g[k]-delta[k],2) < 20.)){
            MultiLinearFunction psi_i_L(D-1), psi_i_U(D-1);
            for(int h = 0; h<psi_i_L.c.size(); h++){
                auto bools = MultiLinearFunction::get_bools(h, D-1);
                bools.insert(bools.begin()+k, 0);
                psi_i_L.c[h] += psi_i[bools]; psi_i_U.c[h] += psi_i[bools];
                bools = MultiLinearFunction::get_bools(h, D-1);
                bools.insert(bools.begin()+k, 1);
                psi_i_L.c[h] += psi_i[bools]*xL[k]; psi_i_U.c[h] += psi_i[bools]*xU[k];
            }
            psitilde.Append(psi_i_L); psitilde.Append(psi_i_U);
            int sign_gk = 0;
            if(g[k] < 0) sign_gk = -1;
            else if(g[k] > 0) sign_gk = 1;
            stilde.Append(sgn_L(sign_gk, s_pruned[i], S)); stilde.Append(sgn_U(sign_gk, s_pruned[i], S));
        }
        else {
            if(subdivlevel >= 10) throw Exception ("Saye asked for Subdivlevel 11!");
            vector<double> d(D);
            for(int j=0; j<D; j++) d[j] = abs(xU[j] - xL[j]);
            int k = distance(d.begin(), max_element(d.begin(), d.end()));
            cout << "Dimension for subdivision: " << k << endl;
            Vec<D> xL1, xL2, xU1, xU2;
            for(int j=0; j<D; j++){
                if(j == k){
                    xL1[j] = xL[j]; xU1[j] = 0.5*(xL[j] + xU[j]);
                    xL2[j] = xU1[j]; xU2[j] = xU[j];
                }
                else{
                    xL1[j] = xL[j]; xL2[j] = xL[j];
                    xU1[j] = xU[j]; xU2[j] = xU[j];
                }
            }
            cout << "xL1: " << xL1 << endl << "xU1: " << xU1 << endl << "xL2: " << xL2 << endl << "xU2: " << xU2 << endl;
            integrate_saye(psi_pruned, s_pruned, xL1, xU1, S, order, result, subdivlevel+1);
            integrate_saye(psi_pruned, s_pruned, xL2, xU2, S, order, result, subdivlevel+1);
            return;
        }
    }
    function<IntegrationRule(Vec<D-1>)> ftilde;
    if(S) ftilde = [&] (Vec<D-1> x) {IntegrationRule ir_r; eval_surface_integrand<D>(psi_pruned[0],k, xL[k], xU[k], x, ir_r); return ir_r; };
    else  ftilde = [&] (Vec<D-1> x) {IntegrationRule ir_r; eval_integrand<D>(psi_pruned, s_pruned, k, xL[k], xU[k], x, order, ir_r); return ir_r; };

    Vec<D-1> xLtilde, xUtilde;
    for(int i=0; i<D; i++){
        if(i < k) { xLtilde[i] = xL[i]; xUtilde[i] = xU[i]; }
        else if (i > k) { xLtilde[i-1] = xL[i]; xUtilde[i-1] = xU[i]; };
    }

    IntegrationRule ir;
    integrate_saye<D-1>(psitilde, stilde, xLtilde, xUtilde, false, order, ir);
    for(auto ip: ir){
        auto ir_l = ftilde(ip.Point());
        for(auto ip2: ir_l) {
            result.Append(IntegrationPoint(ip2.Point(), ip.Weight()*ip2.Weight()));
        }
    }
    //cout << "The Intrule in the integrate_saye function" << endl;
    //cout << result << endl;
  }

  template<>
  void integrate_saye<1>(Array<MultiLinearFunction>& psi, Array<int>& s, Vec<1> xL, Vec<1> xU, bool S, int order, IntegrationRule& result, int subdivlevel) {
    eval_integrand<1>(psi, s, 0, xL[0], xU[0], {}, order, result);
  }

  double DebugSaye(int s_dt, bool IR_mode){
    MultiLinearFunction phi(2);
    //phi[{0,0}] = 1;
    //phi[{0,1}] = -2;
    //phi[{1,0}] = -2;
    //phi[{1,1}] = -1./3;
    phi.FromLsetVals(Vec<4>{1,-1,-3,-1});
    std::function<double(Vec<2>)> f = [](Vec<2> x) -> double {return 1;};

    //cout << phi.get_largest_abs_on_hyperrect<2>({0.,0.}, {2.,2.}) << endl;

    /*cout << "Get_Tensor_Product_IR testing: " << endl;
    IntegrationRule ir;
    Get_Tensor_Product_IR<3>(8, Vec<3>{0,0,0}, Vec<3>{1,1,1}, ir);*/

    //return eval_surface_integrand<2>(phi, 0, 0.,1., Vec<1>{0.25}, f);
    //return eval_integrand<2>({phi}, {-1}, 0, 0.,1., Vec<1>{0.25}, f, 1);

    Array<MultiLinearFunction> phis{phi}; Array<int> sis{s_dt};
    //bool IR_mode = true;
    Array<int> sis2{0};

    if (IR_mode){
        IntegrationRule irIF; integrate_saye(phis ,sis2, Vec<2>{0.,0.}, Vec<2>{1.,1.}, true, 1, irIF);
        cout << "Getting the IF integration rule: " << irIF << endl;
        double I = 0;
        IntegrationRule ir; integrate_saye(phis, sis, Vec<2>{0.,0.}, Vec<2>{1.,1.}, false, 1, ir);
        for(auto ip:ir) {
            I += ip.Weight()*f(ip.Point());
            cout << ip << endl;
        }
        return I;
    }
    else {
        double II = integrate_saye(phis, sis2, Vec<2>{0.,0.}, Vec<2>{1.,1.}, f, true, 1);
        cout << "Interface integration: " << endl;
        cout << II << endl;

        return integrate_saye(phis, sis, Vec<2>{0.,0.}, Vec<2>{1.,1.}, f, false, 1);
    }
  }

  void MultiLinearFunction::FromLsetVals(FlatVector<> lsetvals){
      if(D == 2){
          operator[]({0,0}) = lsetvals[0]; operator[]({1,0}) = lsetvals[1]-operator[]({0,0}), operator[]({0,1}) = lsetvals[3] - operator[]({0,0}), operator[]({1,1}) = lsetvals[2] - operator[]({0,0}) - operator[]({1,0}) - operator[]({0,1});
      }
      else{
          operator[]({0,0,0}) = lsetvals[0];
          operator[]({1,0,0}) = lsetvals[1]-lsetvals[0];
          operator[]({0,1,0}) = lsetvals[3] - lsetvals[0];
          operator[]({0,0,1}) = lsetvals[4] - lsetvals[0];
          operator[]({1,1,0}) = lsetvals[2] - operator[]({1,0,0}) - operator[]({0,1,0}) - operator[]({0,0,0});
          operator[]({1,0,1}) = lsetvals[5] - operator[]({1,0,0}) - operator[]({0,0,1}) - operator[]({0,0,0});
          operator[]({0,1,1}) = lsetvals[7] - operator[]({0,1,0}) - operator[]({0,0,1}) - operator[]({0,0,0});
          operator[]({1,1,1}) = lsetvals[6] - operator[]({1,1,0}) - operator[]({1,0,1}) - operator[]({0,1,1}) - operator[]({1,0,0}) - operator[]({0,1,0}) - operator[]({0,0,1}) - operator[]({0,0,0});
      }
  }

  void SayeCutElementGeometry::GetIntegrationRule(int order, DOMAIN_TYPE dt, IntegrationRule &intrule){
      levelset = MultiLinearFunction(D);
      levelset.FromLsetVals(lset);

      Array<MultiLinearFunction> phis{levelset}; Array<int> sis{0}; bool S = true;
      if(dt == POS) { sis[0] = 1; S = false; }
      else if(dt == NEG) { sis[0] = -1; S = false; }

      if(D == 2){
          integrate_saye<2>(phis, sis, Vec<2>{0.,0.}, Vec<2>{1.,1.}, S, order, intrule);
      }
      else if (D == 3) {
          integrate_saye<3>(phis, sis, Vec<3>{0.,0.,0.}, Vec<3>{1.,1.,1.}, S, order, intrule);
      }
      else {
          throw Exception( "This Dim is not supported yet in SayeCutElementGeometry");
      }
      //cout << "The Intrule: " << endl << intrule << endl;
  }
  Vec<3> SayeCutElementGeometry::GetNormal(const Vec<3>& p) const{
      Vec<3> n;
      if (D == 3){
          n = levelset.get_grad(p);
          n /= L2Norm(n);
      }
      else if(D == 2){
          Vec<2> n_red = levelset.get_grad(Vec<2>{p[0], p[1]});
          n_red /= L2Norm(n_red);
          n[0] = n_red[0]; n[1] = n_red[1]; n[2] = n_red[2];
      }
      else {
          throw Exception("Unsupported Dim in SayeCutElementGeometry::GetNormal");
      }
      return n;
  }

  template<unsigned int D>
  void TransformQuadUntrafoToIRInterface(const IntegrationRule & quad_untrafo, const ElementTransformation & trafo, const CutElementGeometry & geom, IntegrationRule * ir_interface){
      for (int i = 0; i < quad_untrafo.Size(); ++i)
      {
          MappedIntegrationPoint<D,D> mip(quad_untrafo[i],trafo);
          Mat<D,D> Finv = mip.GetJacobianInverse();

          Vec<3> normal = Trans(Finv) * geom.GetNormal(quad_untrafo[i].Point()) ;
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
                                                     LocalHeap & lh, bool use_saye)
  {
    static Timer t ("NewStraightCutIntegrationRule");
    static Timer timercutgeom ("NewStraightCutIntegrationRule::CheckIfCutFast");
    static Timer timermakequadrule("NewStraightCutIntegrationRule::MakeQuadRule");

    RegionTimer reg(t);

    int DIM = trafo.SpaceDim();

    auto et = trafo.GetElementType();

    if ((et != ET_TRIG)&&(et != ET_TET)&&(et != ET_QUAD)&&(et != ET_HEX)){
      cout << "Element Type: " << et << endl;
      throw Exception("only trigs, tets and quads for now");
    }

    timercutgeom.Start();
    auto element_domain = CheckIfStraightCut(cf_lset_at_element);
    timercutgeom.Stop();

    timermakequadrule.Start();
    CutSimplexElementGeometry geom(cf_lset_at_element, et, lh);
    CutQuadElementGeometry geom_quad(cf_lset_at_element, et, lh);
    SayeCutElementGeometry geom_quad_saye(cf_lset_at_element, et, lh);
    IntegrationRule quad_untrafo;

    if (element_domain == IF)
    {
      static Timer timer1("StraightCutElementGeometry::Load+Cut");
      timer1.Start();
      if((et == ET_QUAD)||(et == ET_HEX)) {
          if(!use_saye) geom_quad.GetIntegrationRule(intorder, dt, quad_untrafo);
          else geom_quad_saye.GetIntegrationRule(intorder, dt, quad_untrafo);
      }
      else geom.GetIntegrationRule(intorder, dt, quad_untrafo);
      timer1.Stop();
    }

    const IntegrationRule* ir = nullptr;

    timermakequadrule.Stop();

    if (element_domain == IF) // there is a cut on the current element
    {
      if (dt == IF)
      {
        auto ir_interface  = new (lh) IntegrationRule(quad_untrafo.Size(),lh);
        if (DIM == 2){
            if(et == ET_QUAD){
                if(!use_saye) TransformQuadUntrafoToIRInterface<2>(quad_untrafo, trafo, geom_quad, ir_interface);
                else TransformQuadUntrafoToIRInterface<2>(quad_untrafo, trafo, geom_quad_saye, ir_interface);
            }
            else TransformQuadUntrafoToIRInterface<2>(quad_untrafo, trafo, geom, ir_interface);
        }
        else{
            if(et == ET_HEX){
                if(!use_saye) TransformQuadUntrafoToIRInterface<3>(quad_untrafo, trafo, geom_quad, ir_interface);
                else TransformQuadUntrafoToIRInterface<3>(quad_untrafo, trafo, geom_quad_saye, ir_interface);
            }
            else TransformQuadUntrafoToIRInterface<3>(quad_untrafo, trafo, geom, ir_interface);
        }
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
} // end of namespace
