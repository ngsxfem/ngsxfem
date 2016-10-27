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

  SimpleX StraightCutElementGeometry::Cut(SimpleX s){
    /*cout << "Cut was called for a Simplex with " << s.size() << " Elements: ";
    for(auto t:s) cout << t << "\t";
    cout << endl;*/
    if(s.size() == 2){
        svs.push_back(svs[s[0]]+(lset[s[0]]/(lset[s[0]]-lset[s[1]]))*(svs[s[1]]-svs[s[0]]));
        return {{svs.size()-1}};
    }
    else if(s.size() >= 3){
        set<int> cut_points;
        for(int i_remove =0; i_remove<s.size(); i_remove++){
            SimpleX red = s;
            red.erase(red.begin()+i_remove);
            FlatVector<> lset_red(red.size(), lh);
            for(int i=0; i<red.size(); i++) lset_red[i] = lset[red[i]];
            if (CheckIfStraightCut(lset_red) == IF){
                SimpleX s_i = Cut(red);
                for(int t: s_i) cut_points.insert(t);
                //cut_points.insert(cut_points.begin(), s_i.begin(), s_i.end());
                /*cout << "Inserted the following SimpleX: ";
                for(auto t: s_i) cout << t << " ";
                cout << endl;*/
            }
        }
        /*if(s.size() > 3) RemoveDuplicates(cut_points);
        cout << "The Vector cut_points: " << endl;
        for(auto p : cut_points) cout << p << "\t [ = " << svs[p] << "]" << endl;
        cout << endl << endl;*/
        SimpleX cut_points_vec; cut_points_vec.assign(cut_points.begin(), cut_points.end());
        return cut_points_vec;
    }
  }

  double StraightCutElementGeometry::MeasureSimplVol(SimpleX s){
      if(s.size()==2) return L2Norm(Vec<3>(svs[s[1]]-svs[s[0]]));
      else if(s.size()==3) return L2Norm(Cross(Vec<3>(svs[s[2]]-svs[s[0]]), Vec<3>(svs[s[1]]-svs[s[0]])));
      else if(s.size()==4) return abs(Determinant<3>(Vec<3>(svs[s[3]]-svs[s[0]]), Vec<3>(svs[s[2]]-svs[s[0]]), Vec<3>(svs[s[1]]-svs[s[0]])));
      else throw Exception("Calc the Volume of this type of Simplex not implemented!");
  }

  void StraightCutElementGeometry::AppendIntegrationRuleOnSimpl(SimpleX s, int order, IntegrationRule& intrule){
      double trafofac = MeasureSimplVol(s);

      IntegrationRule ir_ngs;
      if(s.size() == 2) ir_ngs = SelectIntegrationRule(ET_SEGM, order);
      else if(s.size() == 3) ir_ngs = SelectIntegrationRule(ET_TRIG, order);
      else if(s.size() == 4) ir_ngs = SelectIntegrationRule (ET_TET, order);

      IntegrationRule ir_new;
      for (auto ip : ir_ngs) {
        Vec<3> point(0.0); double originweight = 1.0;
        for (int m = 0; m < s.size()-1 ;++m) originweight -= ip(m);
        point = originweight * (svs[s[0]]);
        for (int m = 0; m < s.size()-1 ;++m)
          point += ip(m) * (svs[s[m+1]]);
        double weight = ip.Weight() * trafofac;
        intrule.AddIntegrationPoint(IntegrationPoint(point, weight));
        //cout << point << "\t W: " << weight << endl;
      }
  }

  void StraightCutElementGeometry::LoadBaseSimplexFromElementTopology() {
      const POINT3D * verts = ElementTopology::GetVertices(et);

      for(int i=0; i<ElementTopology::GetNVertices(et); i++)
          svs.push_back({verts[i][0], verts[i][1], verts[i][2]});

      if(et == ET_TRIG) { // || et == ET_TET){ -- To be added after the cutting form ET_TET is finished
          vector<int> BaseSimplex(D+1); iota(BaseSimplex.begin(), BaseSimplex.end(), 0);
          simplices.push_back(BaseSimplex);
      }
      else throw Exception("Error in LoadBaseSimplexFromElementTopology() - ET_TYPE not supported yet!");
  }

  void StraightCutElementGeometry::CalcNormal(SimpleX s_cut){
      if(D == 2) normal = {Vec<3>(svs[s_cut[1]]-svs[s_cut[0]])[1], - Vec<3>(svs[s_cut[1]]-svs[s_cut[0]])[0], 0}; //Cross(Vec<3>(svs[s_cut[1]]-svs[s_cut[0]]), Vec<3>(0,0,1));
      else if(D == 3) normal = Cross(Vec<3>(svs[s_cut[2]] - svs[s_cut[0]]), Vec<3>(svs[s_cut[1]] - svs[s_cut[0]]));
      normal /= L2Norm(normal);
      if ((InnerProduct(Vec<3>(svs[0]-svs[s_cut[0]]),normal) < 0)^(lset[0] > 0))
          normal *= -1.;
      //cout << "normal: " << normal << endl;
  }

  void StraightCutElementGeometry::CutBaseSimplex(DOMAIN_TYPE dt){
      SimpleX s_cut = Cut(simplices[0]);

      simplices.clear();
      if(dt == IF) {
          simplices.push_back(s_cut);
          CalcNormal(s_cut);
      }
      else {
          vector<int> relevant_base_simplex_vertices;
          for(int i=0; i<D+1; i++)
              if( ((dt == POS) &&(lset[i] > 0)) || ((dt == NEG) &&(lset[i] < 0)))
                  relevant_base_simplex_vertices.push_back(i);
          if(relevant_base_simplex_vertices.size() == 1){ //Triangle is cut to a triangle || Tetraeder to a tetraeder
              vector<int> s = s_cut;
              s.push_back(relevant_base_simplex_vertices[0]);
              simplices.push_back(s);
          }
          else if((relevant_base_simplex_vertices.size() == 2) && (D==2)){ //Triangle is cut to a quad
              vector<int> s1, s2; s1 = relevant_base_simplex_vertices; s2 = s_cut;
              s1.push_back(s_cut[1]); s2.push_back(relevant_base_simplex_vertices[1]); //The right indices follow from the cutting order
              simplices.push_back(s1);
              simplices.push_back(s2);
          }
          else {
              throw Exception("Cutting this part of a tetraeder is not implemented yet!");
          }
      }
  }

  void StraightCutElementGeometry::GetIntegrationRule(int order, IntegrationRule& intrule){
      for(auto s : simplices) AppendIntegrationRuleOnSimpl(s, order, intrule);
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

    if (et != ET_TRIG)
      throw Exception("only trigs for now");

    timercutgeom.Start();
    auto element_domain = CheckIfStraightCut(cf_lset_at_element);
    timercutgeom.Stop();

    //cout << "We are at Element Nr. " << trafo.GetElementNr() << endl;

    /*
    if(trafo.GetElementNr() == 907) {
        cout << "Last Element!! Debugging stuff:" << endl << endl;
        FlatVector<> lset(4,lh); lset = Vec<4>{-1,-1,-1,1};
        StraightCutElementGeometry geom(lset, et, lh);
        Simpl s {{0,0,0}, {1,0,0}, {0,1,0},{0,0,1}};
        geom.svs = s; geom.simplices = {{0,1,2,3}};
        cout << "geom.lset: " << geom.lset << endl;
        geom.Cut(geom.simplices[0], lh);
        cout << "End of debugging stuff!" << endl;
    }*/

    timermakequadrule.Start();
    StraightCutElementGeometry geom(cf_lset_at_element, et, lh);

    if (element_domain == IF)
    {
      static Timer timer1("StraightCutElementGeometry::Load+Cut"); timer1.Start();
      geom.LoadBaseSimplexFromElementTopology();
      geom.CutBaseSimplex(dt);
      timer1.Stop();
    }

    IntegrationRule quad_untrafo;
    geom.GetIntegrationRule(intorder, quad_untrafo);
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
              const double weight = interface_quad.weights[i] * len;

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
              const double weight = interface_quad.weights[i] * len;

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
