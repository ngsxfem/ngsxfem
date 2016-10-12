#include "straightcutrule.hpp"

namespace xintegration
{

  // a poor mans try to speed things up (only gave a factor of 2) - can be optimized further...
  DOMAIN_TYPE StraightCutDomain(shared_ptr<CoefficientFunction> cf_lset,
                                const ElementTransformation & trafo,
                                LocalHeap & lh)
  { // determine interface_normal_ref
    static Timer timera ("HasStraightCut");
    RegionTimer reg (timera);

    bool haspos = false;
    bool hasneg = false;

    const POINT3D * verts = ElementTopology::GetVertices(trafo.GetElementType());
    const int nv =  ElementTopology::GetNVertices(trafo.GetElementType());
    for (int v = 0; v < nv; ++v)
    {
      IntegrationPoint vip(&verts[v][0],0);
      MappedIntegrationPoint<2,2> mip(vip,trafo);
      const double tmp = cf_lset->Evaluate(mip);
      if (tmp >= 0.0)
        haspos = true;
      else
        hasneg = true;
    }
    if (hasneg && haspos)
      return IF;
    else if (hasneg)
      return NEG;
    else
      return POS;
  }

  // integration rules that are returned assume that a scaling with mip.GetMeasure() gives the
  // correct weight on the "physical" domain (note that this is not a natural choicefor interface integrals)
  const IntegrationRule * StraightCutIntegrationRule(shared_ptr<CoefficientFunction> cf_lset,
                                                     const ElementTransformation & trafo,
                                                     DOMAIN_TYPE dt,
                                                     int intorder,
                                                     LocalHeap & lh)
  {
    static Timer t ("StraightCutRule");
    static Timer timercutgeom ("StraightCutRule::MakeQuadRule");

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

    auto element_domain = StraightCutDomain(cf_lset,trafo,lh);

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
      // DOMAIN_TYPE element_domain =
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
