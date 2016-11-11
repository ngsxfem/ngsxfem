#include "straightcutrule.hpp"

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

  void StraightCutElementGeometry::LoadBaseSimplexFromElementTopology() {
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

  void StraightCutElementGeometry::CalcNormal(const Polytope &base_simplex){
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

  void StraightCutElementGeometry::CutBaseSimplex(DOMAIN_TYPE dt){
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

  void StraightCutElementGeometry::GetIntegrationRule(int order, DOMAIN_TYPE dt, IntegrationRule &intrule){
      LoadBaseSimplexFromElementTopology();
      CutBaseSimplex(dt);

      int ref_ir_ngs_size; int j =0;
      for(int i=0; i<simplices.Size(); i++) {
        double trafofac = MeasureSimplVol(simplices[i]);

        IntegrationRule ir_ngs;
        if(simplices[i].Size() == 2) ir_ngs = SelectIntegrationRule(ET_SEGM, order);
        else if(simplices[i].Size() == 3) ir_ngs = SelectIntegrationRule(ET_TRIG, order);
        else if(simplices[i].Size() == 4) ir_ngs = SelectIntegrationRule (ET_TET, order);

        //cout << "Size of Simplices nr. " << i << ": " << simplices[i].Size() << endl;

        if(i == 0){
            intrule.SetSize(simplices.Size()*ir_ngs.Size());
            ref_ir_ngs_size = ir_ngs.Size();
        }
        else if (ir_ngs.Size() != ref_ir_ngs_size) throw Exception("Different sizes for ir_ngs are not supported!");

        for (auto ip : ir_ngs) {
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

  void StraightCutQuadElementGeometry::LoadBaseQuadFromElementTopology(FlatVector<> lset) {
      const POINT3D * verts = ElementTopology::GetVertices(et);

      for(int i=0; i<ElementTopology::GetNVertices(et); i++)
          svs_ptr->Append(make_tuple(Vec<3>{verts[i][0], verts[i][1], verts[i][2]}, lset[i]));

      if(et != ET_QUAD) throw Exception("Error in LoadBaseSimplexFromElementTopology() - ET_TYPE not supported yet!");
  }

  void AppendIntRuleOnScaledQuad(const Polytope & quad, IntegrationRule &intrule, int order){
      Mat<2,2> A(0.0);
      A(0,0) = quad.GetPoint(1)[0] - quad.GetPoint(0)[0];
      A(1,0) = quad.GetPoint(1)[1] - quad.GetPoint(0)[1];
      A(0,1) = quad.GetPoint(3)[0] - quad.GetPoint(0)[0];
      A(1,1) = quad.GetPoint(3)[1] - quad.GetPoint(0)[1];
      IntegrationRule ir_ngs = SelectIntegrationRule(ET_QUAD, order);
      double trafofac = abs(Det(A));
      for(auto ip : ir_ngs){
          Vec<3> p(0.); p = quad.GetPoint(0) + A*ip.Point();
          intrule.Append(IntegrationPoint(p, ip.Weight()*trafofac));
      }
  }

  Vec<3> StraightCutQuadElementGeometry::GetNormal(const Vec<3>& p) const{
      Vec<3> geom_normal(0.);
      geom_normal[0] = a+c*p[1];
      geom_normal[1] = b+c*p[0];
      geom_normal /= L2Norm(geom_normal);
      return geom_normal;
  }


  void StraightCutQuadElementGeometry::GetIntegrationRule(FlatVector<> lset, int order, DOMAIN_TYPE dt, IntegrationRule &intrule){
      LoadBaseQuadFromElementTopology(lset);
      //Ansatz: levelset = a*x+b*y+c*x*y+d
      d = lset[0]; a = lset[1]-d, b = lset[3] - d, c = lset[2] - a- b- d;
      function<double(Vec<3>)> levelset = [this] (Vec<3> p) {return a*p[0]+b*p[1]+c*p[0]*p[1]+d;};

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

      Array<Polytope> Cut_quads; Array<Polytope> Volume_quads;
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

      if(dt == IF){
          for(auto poly: Cut_quads){
              double x0 = poly.GetPoint(0)[0], x1 = poly.GetPoint(1)[0];
              IntegrationRule ir_ngs;
              ir_ngs = SelectIntegrationRule(ET_SEGM, order);

              for (auto ip : ir_ngs) {
                Vec<3> point(0.0);
                point[0] = x0+ip.Point()[0]*(x1-x0);
                double u = a*point[0]+d, v = c*point[0]+b; point[1] = -u/v;
                Vec<3> grad_gamma(0.0); grad_gamma[0] = x1-x0;
                grad_gamma[1] = -(x1-x0)*(a*v - c*u)/(pow(v,2));
                intrule.Append(IntegrationPoint(point, ip.Weight() * L2Norm(grad_gamma)));
              }
          }
      }
      else {
          for(auto poly : Volume_quads) AppendIntRuleOnScaledQuad(poly, intrule, order);
          for(auto poly : Cut_quads){
              double x0 = poly.GetPoint(0)[0], x1 = poly.GetPoint(2)[0];
              function<double(double)> y_ast = [this](double x) -> double {return -(a*x+d)/(c*x+b);};
              auto y0 = y_ast, y1 = y_ast;
              if(((dt == POS)&&(poly.GetLset(0) > 1e-12))||((dt == NEG)&&(poly.GetLset(0) < -1e-12)))
                  y0 = [&poly] (double x) -> double {return poly.GetPoint(0)[1];};
              else
                  y1 = [&poly] (double x) ->double {return poly.GetPoint(2)[1];};

              IntegrationRule ir_ngs = SelectIntegrationRule(ET_SEGM, order);
              Vec<2> scale_f; scale_f[0] = x1-x0;
              for(auto p1: ir_ngs){
                  for(auto p2 : ir_ngs){
                      Vec<3> p(0.); p[0] = x0+p1.Point()[0]*scale_f[0];
                      scale_f[1] = y1(p[0]) - y0(p[0]);
                      p[1] = y0(p[0]) + p2.Point()[0]*scale_f[1];
                      intrule.Append(IntegrationPoint(p, p1.Weight()*p2.Weight()*scale_f[0]*scale_f[1]));
                  }
              }
          }
      }
  }

  template<unsigned int D, typename T>
  void TransformQuadUntrafoToIRInterface(const IntegrationRule & quad_untrafo, const ElementTransformation & trafo, const T & geom, IntegrationRule * ir_interface){
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
  const IntegrationRule * StraightCutIntegrationRule(shared_ptr<CoefficientFunction> cf_lset,
                                                     const FlatVector<> & cf_lset_at_element,
                                                     const ElementTransformation & trafo,
                                                     DOMAIN_TYPE dt,
                                                     int intorder,
                                                     LocalHeap & lh)
  {
    static Timer t ("NewStraightCutIntegrationRule");
    static Timer timercutgeom ("NewStraightCutIntegrationRule::CheckIfCutFast");
    static Timer timermakequadrule("NewStraightCutIntegrationRule::MakeQuadRule");

    RegionTimer reg(t);

    int DIM = trafo.SpaceDim();

    auto et = trafo.GetElementType();

    if ((et != ET_TRIG)&&(et != ET_TET)&&(et != ET_QUAD))
      throw Exception("only trigs, tets and quads for now");

    timercutgeom.Start();
    auto element_domain = CheckIfStraightCut(cf_lset_at_element);
    timercutgeom.Stop();

    timermakequadrule.Start();
    StraightCutElementGeometry geom(cf_lset_at_element, et, lh);
    StraightCutQuadElementGeometry geom_quad(et);
    IntegrationRule quad_untrafo;

    if (element_domain == IF)
    {
      static Timer timer1("StraightCutElementGeometry::Load+Cut");
      timer1.Start();
      if(et == ET_QUAD) geom_quad.GetIntegrationRule(cf_lset_at_element, intorder, dt, quad_untrafo);
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
            if(et == ET_QUAD) TransformQuadUntrafoToIRInterface<2>(quad_untrafo, trafo, geom_quad, ir_interface);
            else TransformQuadUntrafoToIRInterface<2>(quad_untrafo, trafo, geom, ir_interface);
        }
        else TransformQuadUntrafoToIRInterface<3>(quad_untrafo, trafo, geom, ir_interface);
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
