#include "straightcutrule.hpp"

namespace xintegration
{

  // a poor mans try to speed things up (only gave a factor of 2) - can be optimized further...
  DOMAIN_TYPE StraightCutDomain(shared_ptr<CoefficientFunction> cf_lset,
                                const ElementTransformation & trafo,
                                LocalHeap & lh, vector<int>& sign_of_lset_at_vertex, vector<int>& my_vert_idx)
  { // determine interface_normal_ref
    static Timer timera ("StraightCutDomain");
    RegionTimer reg (timera);
    static Timer timerb ("StraightCutDomain::GettingVertices");
    static Timer timerc ("StraightCutDomain::GetMip");
    static Timer timerd ("StraightCutDomain::CfEval");

    bool haspos = false;
    bool hasneg = false;

    timerb.Start();
    const POINT3D * verts = ElementTopology::GetVertices(trafo.GetElementType());
    const int nv =  ElementTopology::GetNVertices(trafo.GetElementType());
    timerb.Stop();

    for (int v = 0; v < nv; ++v)
    {
      if(sign_of_lset_at_vertex[my_vert_idx[v]] == -10){
        timerc.Start();
        IntegrationPoint vip(&verts[v][0],0);
        MappedIntegrationPoint<2,2> mip(vip,trafo);
        timerc.Stop();
        timerd.Start();
        const double tmp = cf_lset->Evaluate(mip);//cf_lset->Evaluate()
        timerd.Stop();
        sign_of_lset_at_vertex[my_vert_idx[v]] = (tmp >= 0.0);
        if (tmp >= 0.0)
            haspos = true;
        else
            hasneg = true;
      }
      else {
          if(sign_of_lset_at_vertex[my_vert_idx[v]] == 1) haspos = true;
          else if(sign_of_lset_at_vertex[my_vert_idx[v]] == 0) hasneg = true;
          else cout << "!! - Encountered bad sign_of_lset_at_vertey-Flag!" << endl;
      }
      if(haspos && hasneg) break;
    }
    if (hasneg && haspos)
      return IF;
    else if (hasneg)
      return NEG;
    else
      return POS;
  }

  // integration rules that are returned assume that a scaling with mip.GetMeasure() gives the
  // correct weight on the "physical" domain (note that this is not a natural choice for interface integrals)
  const IntegrationRule * StraightCutIntegrationRule(shared_ptr<CoefficientFunction> cf_lset,
                                                     const ElementTransformation & trafo,
                                                     DOMAIN_TYPE dt,
                                                     int intorder,
                                                     LocalHeap & lh, vector<int>& sign_of_lset_at_vertex,
                                                     vector<int>& my_vert_idx)
  {
    static Timer t ("StraightCutIntegrationRule");
    static Timer timercutgeom ("StraightCutIntegrationRule::MakeQuadRule");

    int subdivlvl = 0;

    RegionTimer reg(t);

    int DIM = trafo.SpaceDim();
    auto lset_eval
      = ScalarFieldEvaluator::Create(DIM,*cf_lset,trafo,lh);
      // tstart.Stop();
    timercutgeom.Start();

    auto et = trafo.GetElementType();

    if (et != ET_TRIG)
      throw Exception("only trigs for now");

    shared_ptr<XLocalGeometryInformation> xgeom = nullptr;
    CompositeQuadratureRule<2> cquad2d;
    CompositeQuadratureRule<3> cquad3d;


    auto element_domain = StraightCutDomain(cf_lset,trafo,lh,sign_of_lset_at_vertex, my_vert_idx);

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
      xgeom->MakeQuadRule();
    }
    timercutgeom.Stop();


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
