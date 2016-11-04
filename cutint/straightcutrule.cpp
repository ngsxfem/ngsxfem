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

  Polytope StraightCutElementGeometry::CalcCutPointLineUsingLset(const Polytope &s){
      if((s.D != 1) ||(s.Size() != 2)) throw Exception("You called the cut-a-line function with a Polytope which is not a line!");

      Vec<3> p = svs[s[0]]+(lset[s[0]]/(lset[s[0]]-lset[s[1]]))*(svs[s[1]]-svs[s[0]]);
      svs.Append(p);
      return Polytope({svs.Size()-1}, 0);
  }

  Polytope StraightCutElementGeometry::CalcCutPolytopeUsingLset(const Polytope &s){
      Array<int> cut_points;
      if(((s.Size() == 3)&&(s.D == 2))||((s.Size() == 4)&&(s.D == 3))){
        for(int i:s) {
            for(int j:s){
                if(i < j) {
                    if(lset[i]*lset[j] < -1e-10){
                        Polytope p = CalcCutPointLineUsingLset(Polytope({i, j}, 1));
                        cut_points.Append(p[0]);
                    }
                }
            }
        }
    }
    else {
        throw Exception("You tried to cut a Polytope which is not a simplex.");
        return Polytope();
    }
    return Polytope(cut_points, s.D-1);
  }

  double StraightCutElementGeometry::MeasureSimplVol(const Polytope &s){
      if(s.Size()==2) return L2Norm(Vec<3>(svs[s[1]]-svs[s[0]]));
      else if(s.Size()==3) return L2Norm(Cross(Vec<3>(svs[s[2]]-svs[s[0]]), Vec<3>(svs[s[1]]-svs[s[0]])));
      else if(s.Size()==4) return abs(Determinant<3>(Vec<3>(svs[s[3]]-svs[s[0]]), Vec<3>(svs[s[2]]-svs[s[0]]), Vec<3>(svs[s[1]]-svs[s[0]])));
      else throw Exception("Calc the Volume of this type of Simplex not implemented!");
  }

  void StraightCutElementGeometry::LoadBaseSimplexFromElementTopology() {
      const POINT3D * verts = ElementTopology::GetVertices(et);

      for(int i=0; i<ElementTopology::GetNVertices(et); i++)
          svs.Append({verts[i][0], verts[i][1], verts[i][2]});

      if((et == ET_TRIG) || (et == ET_TET)){
          Array<int> BaseSimplex(D+1);
          for(int i=0; i<BaseSimplex.Size(); i++) BaseSimplex[i] = i;
          simplices.Append(Polytope(BaseSimplex, D));
      }
      else throw Exception("Error in LoadBaseSimplexFromElementTopology() - ET_TYPE not supported yet!");
  }

  /*
  void StraightCutElementGeometry::CalcNormal(const Polytope &s_cut){
      if(D == 2) normal = {Vec<3>(svs[s_cut[1]]-svs[s_cut[0]])[1], - Vec<3>(svs[s_cut[1]]-svs[s_cut[0]])[0], 0}; //Cross(Vec<3>(svs[s_cut[1]]-svs[s_cut[0]]), Vec<3>(0,0,1));
      else if(D == 3) normal = Cross(Vec<3>(svs[s_cut[2]] - svs[s_cut[0]]), Vec<3>(svs[s_cut[1]] - svs[s_cut[0]]));
      normal /= L2Norm(normal);
      if ((InnerProduct(Vec<3>(svs[0]-svs[s_cut[0]]),normal) < 0)^(lset[0] > 0))
          normal *= -1.;
      //cout << "normal: " << normal << endl;
  }*/

  void StraightCutElementGeometry::CalcNormal(){
      Vec<3> delta_vec;
      double delta_f;
      Vec<3> grad_f; grad_f = 0;
      for(int i=0; i<D; i++) {
          delta_vec = svs[i]-svs[D];
          delta_f = lset[i]-lset[D];

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

      simplices.DeleteAll();
      if(dt == IF) {
          if(s_cut.Size() == D) simplices.Append(s_cut);
          else if((s_cut.Size() == 4)&&(D==3)){
              simplices.Append(Polytope({s_cut[0],s_cut[1],s_cut[3]},2));
              simplices.Append(Polytope({s_cut[0],s_cut[2],s_cut[3]},2));
          }
          else {
              cout << "s_cut: " << s_cut.ia << endl;
              throw Exception("Bad length of s_cut!");
          }
          CalcNormal();
      }
      else {
          Polytope relevant_base_simplex_vertices;
          for(int i=0; i<D+1; i++)
              if( ((dt == POS) &&(lset[i] > 1e-10)) || ((dt == NEG) &&(lset[i] < -1e-10)))
                  relevant_base_simplex_vertices.Append(i);
          if((relevant_base_simplex_vertices.Size() == 1)){ //Triangle is cut to a triangle || Tetraeder to a tetraeder
              Polytope s(s_cut);
              s.Append(relevant_base_simplex_vertices[0]);
              simplices.Append(s);
          }
          else if((relevant_base_simplex_vertices.Size() == 2) && (D==2)){ //Triangle is cut to a quad
              Polytope s1, s2; s1 = relevant_base_simplex_vertices; s2 = s_cut;
              s1.Append(s_cut[1]); s2.Append(relevant_base_simplex_vertices[0]); //The right indices follow from the cutting order
              simplices.Append(s1);
              simplices.Append(s2);
          }
          else if((relevant_base_simplex_vertices.Size() == 2) && (D==3)) { //Tetraeder is cut to several tetraeder
              Polytope s1, s2, s3; s1 = s_cut; s1.ia[0] = relevant_base_simplex_vertices[1];
              s2 = relevant_base_simplex_vertices; s2.Append(s_cut[1]); s2.Append(s_cut[2]);
              s3 = s_cut; s3.ia[3] = relevant_base_simplex_vertices[0];
              simplices.Append(s1);
              simplices.Append(s2);
              simplices.Append(s3);
          }
          else if((relevant_base_simplex_vertices.Size() == 3) && (D == 3)){
              Polytope s1, s2, s3; s1 = s_cut; s1.Append(relevant_base_simplex_vertices[2]);
              s2 = relevant_base_simplex_vertices; s2.Append(s_cut[1]);
              s3 = s_cut; s3.ia[2] = relevant_base_simplex_vertices[0]; s3.Append(relevant_base_simplex_vertices[2]);
              simplices.Append(s1);
              simplices.Append(s2);
              simplices.Append(s3);
          }
          else {
              throw Exception("Cutting this part of a tetraeder is not implemented yet!");
          }
      }
  }

  void StraightCutElementGeometry::GetIntegrationRule(int order, IntegrationRule &intrule){
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
            point = originweight * (svs[simplices[i][0]]);
          for (int m = 0; m < simplices[i].Size()-1 ;++m)
            point += ip(m) * (svs[simplices[i][m+1]]);
          intrule[j] = IntegrationPoint(point, ip.Weight() * trafofac);
          j++;
        }
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

    if ((et != ET_TRIG)&&(et != ET_TET))
      throw Exception("only trigs and tets for now");

    timercutgeom.Start();
    auto element_domain = CheckIfStraightCut(cf_lset_at_element);
    timercutgeom.Stop();

    //cout << "We are at Element Nr. " << trafo.GetElementNr() << endl;

    /*
    if(trafo.GetElementNr() == 1371) {
        cout << "Last Element!! Debugging stuff:" << endl << endl;
        FlatVector<> lset(4,lh); lset = Vec<4>{-1,-1,1,1};
        StraightCutElementGeometry geom(lset, ET_TET, lh);
        geom.svs = {{1,0,0}, {0,1,0}, {0,0,1}, {0,0,0}}; geom.simplices = {Polytope({0,1,2,3},3)};
        cout << "geom.lset: " << geom.lset << endl;
        //geom.CalcCutPolytopeUsingLset(geom.simplices[0]);
        geom.CutBaseSimplex(IF);

        cout << "End of debugging stuff!" << endl;
    }*/

    timermakequadrule.Start();
    StraightCutElementGeometry geom(cf_lset_at_element, et, lh);
    IntegrationRule quad_untrafo;

    if (element_domain == IF)
    {
      static Timer timer1("StraightCutElementGeometry::Load+Cut"); timer1.Start();
      geom.LoadBaseSimplexFromElementTopology();
      geom.CutBaseSimplex(dt);
      geom.GetIntegrationRule(intorder, quad_untrafo);
      timer1.Stop();
    }

    const IntegrationRule* ir = nullptr;

    timermakequadrule.Stop();

    if (element_domain == IF) // there is a cut on the current element
    {
      if (dt == IF)
      {
        if (DIM == 2)
        {
          IntegrationRule * ir_interface  = new (lh) IntegrationRule(quad_untrafo.Size(),lh);
          for (int i = 0; i < quad_untrafo.Size(); ++i)
          {
            MappedIntegrationPoint<2,2> mip(quad_untrafo[i],trafo);
            Mat<2,2> Finv = mip.GetJacobianInverse();

            Vec<3> normal = Trans(Finv) * geom.normal ;
            const double weight = quad_untrafo[i].Weight() * L2Norm(normal);

            (*ir_interface)[i] = IntegrationPoint (quad_untrafo[i].Point(), weight);
          }
          ir = ir_interface;
        }
        else
        {
          IntegrationRule * ir_interface  = new (lh) IntegrationRule(quad_untrafo.Size(),lh);
          for (int i = 0; i < quad_untrafo.Size(); ++i)
          {
            MappedIntegrationPoint<3,3> mip(quad_untrafo[i],trafo);
            Mat<3,3> Finv = mip.GetJacobianInverse();

            Vec<3> normal = Trans(Finv) * geom.normal;
            const double weight = quad_untrafo[i].Weight() * L2Norm(normal);

            (*ir_interface)[i] = IntegrationPoint (quad_untrafo[i].Point(), weight);
          }
          ir = ir_interface;
        }
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

  const IntegrationRule * StraightCutIntegrationRuleOld(shared_ptr<CoefficientFunction> cf_lset,
                                                       const FlatVector<> & cf_lset_at_element,
                                                       const ElementTransformation & trafo,
                                                       DOMAIN_TYPE dt,
                                                       int intorder,
                                                       LocalHeap & lh)
    {
      static Timer t ("StraightCutIntegrationRule");
      static Timer timercutgeom ("StraightCutIntegrationRule::CheckIfCutFast");
      static Timer timermakequadrule("StraightCutIntegrationRule::MakeQuadRule");

      //cout << "OLD: We are at Element Nr. " << trafo.GetElementNr() << endl;

      int subdivlvl = 0;

      RegionTimer reg(t);

      int DIM = trafo.SpaceDim();
      auto lset_eval
        = ScalarFieldEvaluator::Create(DIM,*cf_lset,trafo,lh);
      timercutgeom.Start();

      auto et = trafo.GetElementType();

      if (et != ET_TRIG)
        throw Exception("only trigs for now");

      shared_ptr<XLocalGeometryInformation> xgeom = nullptr;
      CompositeQuadratureRule<2> cquad2d;
      CompositeQuadratureRule<3> cquad3d;

      //auto element_domain = StraightCutDomain(cf_lset,trafo,lh);
      auto element_domain = CheckIfStraightCut(cf_lset_at_element);
      timercutgeom.Stop();

      if (element_domain == IF)
      {
        if (DIM == 2)
          xgeom = XLocalGeometryInformation::Create(et, ET_POINT,
                                                    *lset_eval, cquad2d, lh,
                                                    intorder, 0, subdivlvl, 0);
        else
          xgeom = XLocalGeometryInformation::Create(et, ET_POINT,
                                                    *lset_eval, cquad3d, lh,
                                                     intorder, 0, subdivlvl, 0);
        xgeom->cf_lset_at_element.AssignMemory(cf_lset_at_element.Size(), cf_lset_at_element.Data());

        timermakequadrule.Start();
        xgeom->MakeQuadRule();
        //xgeom->MakeQuadRuleFast();
        timermakequadrule.Stop();
      }

      const IntegrationRule* ir = nullptr;

      if (element_domain == IF) // there is a cut on the current element
      {
        if (dt == IF)
        {
          if (DIM == 2)
          {
            const QuadratureRuleCoDim1<2> & interface_quad(cquad2d.GetInterfaceRule());
            IntegrationRule * ir_interface  = new (lh) IntegrationRule(interface_quad.Size(),lh);
            for (int i = 0; i < interface_quad.Size(); ++i)
            {
              IntegrationPoint ip(&interface_quad.points[i](0),interface_quad.weights[i]);
              MappedIntegrationPoint<2,2> mip(ip,trafo);

              Mat<2,2> Finv = mip.GetJacobianInverse();
              const double absdet = mip.GetMeasure();

              Vec<2> nref = interface_quad.normals[i];
              Vec<2> normal = absdet * Trans(Finv) * nref ;
              double len = L2Norm(normal);

              (*ir_interface)[i] = IntegrationPoint (&interface_quad.points[i](0),interface_quad.weights[i] * len / mip.GetMeasure());
              ir = ir_interface;
            }
          }
          else
          {
            const QuadratureRuleCoDim1<3> & interface_quad(cquad3d.GetInterfaceRule());
            IntegrationRule * ir_interface  = new (lh) IntegrationRule(interface_quad.Size(),lh);
            for (int i = 0; i < interface_quad.Size(); ++i)
            {
              IntegrationPoint ip(&interface_quad.points[i](0),interface_quad.weights[i]);
              MappedIntegrationPoint<3,3> mip(ip,trafo);

              Mat<3,3> Finv = mip.GetJacobianInverse();
              const double absdet = mip.GetMeasure();

              Vec<3> nref = interface_quad.normals[i];
              Vec<3> normal = absdet * Trans(Finv) * nref ;
              double len = L2Norm(normal);

              (*ir_interface)[i] = IntegrationPoint (&interface_quad.points[i](0),interface_quad.weights[i] * len / mip.GetMeasure());
              ir = ir_interface;
            }
          }
        }
        else
        {
          if (DIM == 2)
          {
            const QuadratureRule<2> & domain_quad = cquad2d.GetRule(dt);
            auto ir_domain = new (lh) IntegrationRule (domain_quad.Size(),lh);
            for (int i = 0; i < ir_domain->Size(); ++i)
              (*ir_domain)[i] = IntegrationPoint (&domain_quad.points[i](0),domain_quad.weights[i]);
            ir = ir_domain;
          }
          else
          {
            const QuadratureRule<3> & domain_quad = cquad3d.GetRule(dt);
            auto ir_domain = new (lh) IntegrationRule (domain_quad.Size(),lh);
            for (int i = 0; i < ir_domain->Size(); ++i)
              (*ir_domain)[i] = IntegrationPoint (&domain_quad.points[i](0),domain_quad.weights[i]);
            ir = ir_domain;
          }
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
