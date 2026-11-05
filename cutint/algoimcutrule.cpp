#include "algoimcutrule.hpp"
#include "algoim_wrapper.hpp"

namespace xintegration
{
  static_assert(NEG == 0 && POS == 1 && IF == 2,
                "algoim_wrapper::CutQuadrature expects NEG = 0, POS = 1, IF = 2");

  namespace
  {
    template <int D>
    tuple<const IntegrationRule *, Array<double>>
    T_AlgoimCutIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                                const ElementTransformation & trafo,
                                DOMAIN_TYPE dt,
                                int intorder,
                                int degree,
                                int subdivlvl,
                                LocalHeap & lh)
    {
      const ELEMENT_TYPE et = trafo.GetElementType();
      const bool simplex = (et == ET_TRIG) || (et == ET_TET);

      // level set values in the interpolation points in [0,1]^D (for simplices,
      // these points are also outside of the reference element, where the
      // level set function is evaluated by polynomial extension)
      std::vector<double> xi;
      algoim_wrapper::InterpolationPoints(D, degree, xi);
      const int nnodes = xi.size() / D;
      IntegrationRule ir_nodes(nnodes, lh);
      for (int i = 0; i < nnodes; i++)
      {
        Vec<3> p = 0.0;
        for (int d = 0; d < D; d++)
          p(d) = xi[i*D+d];
        ir_nodes[i] = IntegrationPoint(p, 0.0);
      }
      auto & mir_nodes = trafo(ir_nodes, lh);
      FlatMatrix<> lsetvals(nnodes, 1, lh);
      cf_lset->Evaluate(mir_nodes, lsetvals);

      // Algoim's quadrature on the reference element
      const int q = max(1, intorder / 2 + 1);
      std::vector<double> pts, wts, grads;
      auto status = algoim_wrapper::CutQuadrature(D, degree, &lsetvals(0,0), simplex, int(dt), q,
                                                  subdivlvl, pts, wts, grads);

      if (status != algoim_wrapper::CUT)
      {
        const DOMAIN_TYPE element_domain = (status == algoim_wrapper::UNCUT_NEGATIVE) ? NEG : POS;
        if (element_domain != dt)
          return make_tuple(nullptr, Array<double>());
        const IntegrationRule & ir = SelectIntegrationRule(et, intorder);
        Array<double> wei_arr(ir.Size());
        for (int i = 0; i < ir.Size(); i++)
          wei_arr[i] = ir[i].Weight();
        return make_tuple(&ir, std::move(wei_arr));
      }

      const int npts = wts.size();
      if (npts == 0)
        return make_tuple(nullptr, Array<double>());
      auto ir = new (lh) IntegrationRule(npts, lh);
      Array<double> wei_arr(npts);
      for (int i = 0; i < npts; i++)
      {
        Vec<D> x;
        for (int d = 0; d < D; d++)
          x(d) = pts[i*D+d];
        double w = wts[i];
        if (dt == IF)
        {
          // weight w.r.t. the reference surface measure times |F^{-T} n_ref|
          // (F: Jacobian of trafo, n_ref = g / |g|, g the gradient of the level
          // set polynomial), cf. TransformQuadUntrafoToIRInterface
          Vec<D> g;
          for (int d = 0; d < D; d++)
            g(d) = grads[i*D+d];
          const IntegrationPoint ip(x, 0.0);
          MappedIntegrationPoint<D,D> mip(ip, trafo);
          const double gnorm = L2Norm(g);
          w = (gnorm > 0) ? w * L2Norm(Trans(mip.GetJacobianInverse()) * g) / gnorm : 0.0;
        }
        (*ir)[i] = IntegrationPoint(x, w);
        wei_arr[i] = w;
      }
      return make_tuple(ir, std::move(wei_arr));
    }
  }

  bool AlgoimCutRuleAvailable (const ElementTransformation & trafo, DOMAIN_TYPE dt)
  {
    const ELEMENT_TYPE et = trafo.GetElementType();
    if (!((et == ET_SEGM) || (et == ET_TRIG) || (et == ET_QUAD) || (et == ET_TET) || (et == ET_HEX)))
      return false;
    if (dt == IF)
      return (trafo.VB() == VOL) && (et != ET_SEGM)
        && (trafo.SpaceDim() == ElementTopology::GetSpaceDim(et));
    return true;
  }

  int AlgoimLevelsetDegree (const shared_ptr<CoefficientFunction> & cf_lset,
                            const ElementTransformation & trafo,
                            int cf_degree,
                            LocalHeap & lh)
  {
    auto gf = dynamic_pointer_cast<GridFunction>(cf_lset);
    if (gf && (gf->Dimension() == 1))
    {
      const FiniteElement & fel = gf->GetFESpace()->GetFE(trafo.GetElementId(), lh);
      return max(1, fel.Order());
    }
    return max(1, cf_degree);
  }

  tuple<const IntegrationRule *, Array<double>>
  AlgoimCutIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                            const ElementTransformation & trafo,
                            DOMAIN_TYPE dt,
                            int intorder,
                            int cf_degree,
                            int subdivlvl,
                            LocalHeap & lh)
  {
    static Timer timer("AlgoimCutIntegrationRule");
    RegionTimer reg (timer);

    if (!AlgoimCutRuleAvailable(trafo, dt))
      throw Exception("AlgoimCutIntegrationRule: element type / domain type combination not supported");

    const int degree = AlgoimLevelsetDegree(cf_lset, trafo, cf_degree, lh);
    switch (ElementTopology::GetSpaceDim(trafo.GetElementType()))
    {
    case 1: return T_AlgoimCutIntegrationRule<1>(cf_lset, trafo, dt, intorder, degree, subdivlvl, lh);
    case 2: return T_AlgoimCutIntegrationRule<2>(cf_lset, trafo, dt, intorder, degree, subdivlvl, lh);
    case 3: return T_AlgoimCutIntegrationRule<3>(cf_lset, trafo, dt, intorder, degree, subdivlvl, lh);
    default:
      throw Exception("AlgoimCutIntegrationRule: dimension not supported");
    }
  }
}
