#include "straightcutrule.hpp"

namespace xintegration
{
  DOMAIN_TYPE CheckIfStraightCut (FlatVector<> cf_lset_at_element) {
    bool haspos = false;
    bool hasneg = false;

    for (auto v : cf_lset_at_element) {
        if (!haspos && (v >= 0)) haspos = true;
        if (!hasneg && (v <= 0)) hasneg = true;
        if(haspos && hasneg) break;
    }

    if (hasneg && haspos) return IF;
    else if (hasneg) return NEG;
    else return POS;
  }

  SimpleX StraightCutElementGeometry::Cut(const SimpleX &s){
    /*cout << "Cut was called for a Simplex with " << s.size() << " Elements: ";
    for(auto t:s) cout << t << "\t";
    cout << endl;*/
    if(s.Size() == 2){
        Vec<3> p = svs[s[0]]+(lset[s[0]]/(lset[s[0]]-lset[s[1]]))*(svs[s[1]]-svs[s[0]]);
        bool p_in_svs = false; int idx_p_in_svs = 0;
        for(int i=0; i<svs.Size(); i++) if(L2Norm(p- svs[i])<1e-10) {p_in_svs = true; idx_p_in_svs = i; }
        if(p_in_svs) return {idx_p_in_svs};
        else{
            svs.Append(p);
            return {svs.Size()-1};
        }
    }
    else if(s.Size() >= 3){
        Array<int> cut_points;
        for(int i_remove =0; i_remove<s.Size(); i_remove++){
            SimpleX red(s);
            red.DeleteElement(i_remove);
            FlatVector<> lset_red(red.Size(), lh);
            for(int i=0; i<red.Size(); i++) lset_red[i] = lset[red[i]];
            if (CheckIfStraightCut(lset_red) == IF){
                SimpleX s_i = Cut(red);
                for(int t: s_i)
                    if(!cut_points.Contains(t)) cut_points.Append(t);
                //cut_points.insert(cut_points.begin(), s_i.begin(), s_i.end());
                /*cout << "Inserted the following SimpleX: ";
                for(auto t: s_i) cout << t << " ";
                cout << endl;*/
            }
        }
        //if(s.size() > 3) RemoveDuplicates(cut_points);
        /*cout << "The Vector cut_points: " << endl;
        for(auto p : cut_points) cout << p << "\t [ = " << svs[p] << "]" << endl;
        cout << endl << endl;*/
        return cut_points;
    }
  }

  double StraightCutElementGeometry::MeasureSimplVol(const SimpleX &s){
      if(s.Size()==2) return L2Norm(Vec<3>(svs[s[1]]-svs[s[0]]));
      else if(s.Size()==3) return L2Norm(Cross(Vec<3>(svs[s[2]]-svs[s[0]]), Vec<3>(svs[s[1]]-svs[s[0]])));
      else if(s.Size()==4) return abs(Determinant<3>(Vec<3>(svs[s[3]]-svs[s[0]]), Vec<3>(svs[s[2]]-svs[s[0]]), Vec<3>(svs[s[1]]-svs[s[0]])));
      else throw Exception("Calc the Volume of this type of Simplex not implemented!");
  }

  void StraightCutElementGeometry::LoadBaseSimplexFromElementTopology() {
      const POINT3D * verts = ElementTopology::GetVertices(et);

      for(int i=0; i<ElementTopology::GetNVertices(et); i++)
          svs.Append({verts[i][0], verts[i][1], verts[i][2]});

      if((et == ET_TRIG) || (et == ET_TET)){// -- To be added after the cutting form ET_TET is finished
          Array<int> BaseSimplex(D+1);
          for(int i=0; i<BaseSimplex.Size(); i++) BaseSimplex[i] = i;
          simplices.Append(BaseSimplex);
      }
      else throw Exception("Error in LoadBaseSimplexFromElementTopology() - ET_TYPE not supported yet!");
  }

  /*
  void StraightCutElementGeometry::CalcNormal(const SimpleX &s_cut){
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
          //cout << "delta_vec: " << delta_vec << endl;

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
      SimpleX s_cut = Cut(simplices[0]);

      simplices.DeleteAll();
      if(dt == IF) {
          if(s_cut.Size() == D) simplices.Append(s_cut);
          else if((s_cut.Size() == D+1)&&(D==3)){
              simplices.Append({s_cut[1],s_cut[2],s_cut[0]});
              simplices.Append({s_cut[1],s_cut[0],s_cut[3]});
          }
          CalcNormal();
      }
      else {
          Array<int> relevant_base_simplex_vertices;
          for(int i=0; i<D+1; i++)
              if( ((dt == POS) &&(lset[i] > 0)) || ((dt == NEG) &&(lset[i] < 0)))
                  relevant_base_simplex_vertices.Append(i);
          if((relevant_base_simplex_vertices.Size() == 1)&&(D==2)){ //Triangle is cut to a triangle || Tetraeder to a tetraeder
              SimpleX s(s_cut);
              s.Append(relevant_base_simplex_vertices[0]);
              simplices.Append(s);
          }
          else if((relevant_base_simplex_vertices.Size() == 2) && (D==2)){ //Triangle is cut to a quad
              Array<int> s1, s2; s1 = relevant_base_simplex_vertices; s2 = s_cut;
              s1.Append(s_cut[1]); s2.Append(relevant_base_simplex_vertices[1]); //The right indices follow from the cutting order
              simplices.Append(s1);
              simplices.Append(s2);
          }
          /*
          else if((relevant_base_simplex_vertices.Size() == 2) && (D==3)) { //Tetraeder is cut to several tetraeder

          }*/
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

    //if (et != ET_TRIG)
    //  throw Exception("only trigs for now");

    timercutgeom.Start();
    auto element_domain = CheckIfStraightCut(cf_lset_at_element);
    timercutgeom.Stop();

    //cout << "We are at Element Nr. " << trafo.GetElementNr() << endl;

    if(trafo.GetElementNr() == 1371) {
        cout << "Last Element!! Debugging stuff:" << endl << endl;
        FlatVector<> lset(4,lh); lset = Vec<4>{-1,-1,-1,1};
        StraightCutElementGeometry geom(lset, et, lh);
        geom.svs = {{0,0,0}, {1,0,0}, {0,1,0},{0,0,1}}; geom.simplices = {{0,1,2,3}};
        cout << "geom.lset: " << geom.lset << endl;
        geom.Cut(geom.simplices[0]);
        cout << "End of debugging stuff!" << endl;
    }

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
            const double absdet = mip.GetMeasure();

            Vec<2> nref = geom.normal;
            Vec<2> normal = absdet * Trans(Finv) * nref ;
            double len = L2Norm(normal);
            const double weight = quad_untrafo[i].Weight() * len;

            (*ir_interface)[i] = IntegrationPoint (quad_untrafo[i].Point(), weight/ mip.GetMeasure());
            ir = ir_interface;
          }
        }
        else
        {
          IntegrationRule * ir_interface  = new (lh) IntegrationRule(quad_untrafo.Size(),lh);
          for (int i = 0; i < quad_untrafo.Size(); ++i)
          {
            MappedIntegrationPoint<3,3> mip(quad_untrafo[i],trafo);

            Mat<3,3> Finv = mip.GetJacobianInverse();
            const double absdet = mip.GetMeasure();

            Vec<3> nref = geom.normal;
            Vec<3> normal = absdet * Trans(Finv) * nref ;
            double len = L2Norm(normal);
            const double weight = quad_untrafo[i].Weight() * len;

            (*ir_interface)[i] = IntegrationPoint (quad_untrafo[i].Point(),weight * len / mip.GetMeasure());
            ir = ir_interface;
          }
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
