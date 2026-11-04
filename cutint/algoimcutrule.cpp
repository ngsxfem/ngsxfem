#include "algoimcutrule.hpp"
#include "algoim_wrapper.hpp"

namespace xintegration
{
  static_assert(NEG == 0 && POS == 1 && IF == 2,
                "algoim_wrapper::CutQuadrature expects NEG = 0, POS = 1, IF = 2");

  namespace
  {
    // Factor that turns the weight of an interface quadrature point x w.r.t. the
    // reference measure of {phi = 0} (g: gradient of phi at x in reference
    // coordinates) into the weight convention of the cut rules (scaling with the
    // measure of the element transformation gives the physical weight):
    //  * volume elements (S == D): |F^{-T} g| / |g|, F the Jacobian, cf.
    //    TransformQuadUntrafoToIRInterface (for D == 1: points, factor 1/|F|),
    //  * boundary elements (S == D+1), where the interface is of co-dimension 2
    //    in physical space: points (D == 1) have weight 1, curves (D == 2) the
    //    physical length element |F t| (t the unit tangent in reference
    //    coordinates); both divided by the measure of the transformation.
    template <int D, int S>
    double InterfaceWeightFactor (const IntegrationPoint & ip, const ElementTransformation & trafo,
                                  const Vec<D> & g)
    {
      const double gnorm = L2Norm(g);
      MappedIntegrationPoint<D,S> mip(ip, trafo);
      if constexpr (S == D)
        return (gnorm > 0) ? L2Norm(Trans(mip.GetJacobianInverse()) * g) / gnorm : 0.0;
      else if constexpr (D == 1)
        return 1.0 / mip.GetMeasure();
      else
      {
        if (gnorm == 0)
          return 0.0;
        Vec<2> t(g(1) / gnorm, -g(0) / gnorm);
        return L2Norm(mip.GetJacobian() * t) / mip.GetMeasure();
      }
    }

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
      // Algoim uses nested one-dimensional Gauss rules (q points) w.r.t. height functions:
      // * for a polynomial of degree intorder on a domain bounded by planes, the integrand of
      //   the outermost integral has degree intorder + D - 1, hence 2q-1 >= intorder + D - 1,
      // * the (curved) geometry has to be resolved also for low order integrands (e.g. CutInfo
      //   computes cut ratios with intorder 0), hence q > degree.
      const int q = max(degree + 1, (intorder + D + 1) / 2);
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
          Vec<D> g;
          for (int d = 0; d < D; d++)
            g(d) = grads[i*D+d];
          const IntegrationPoint ip(x, 0.0);
          if (trafo.SpaceDim() == D)
            w *= InterfaceWeightFactor<D,D>(ip, trafo, g);
          else if constexpr (D < 3)
            w *= InterfaceWeightFactor<D,D+1>(ip, trafo, g);
        }
        (*ir)[i] = IntegrationPoint(x, w);
        wei_arr[i] = w;
      }
      return make_tuple(ir, std::move(wei_arr));
    }

    template <int D>
    const IntegrationRule *
    T_AlgoimCutFacetIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                                     const ElementTransformation & trafo,
                                     Facet2ElementTrafo & transform, int facetnr,
                                     DOMAIN_TYPE dt, int intorder, int degree, int subdivlvl,
                                     LocalHeap & lh)
    {
      const ELEMENT_TYPE etfacet = transform.FacetType(facetnr);
      const bool simplex = (etfacet == ET_TRIG);

      // level set values in the interpolation points (in facet coordinates)
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
      IntegrationRule & ir_nodes_vol = transform(facetnr, ir_nodes, lh);
      auto & mir_nodes = trafo(ir_nodes_vol, lh);
      FlatMatrix<> lsetvals(nnodes, 1, lh);
      cf_lset->Evaluate(mir_nodes, lsetvals);

      const int q = max(degree + 1, (intorder + D + 1) / 2);
      std::vector<double> pts, wts, grads;
      auto status = algoim_wrapper::CutQuadrature(D, degree, &lsetvals(0,0), simplex, int(dt), q,
                                                  subdivlvl, pts, wts, grads);
      if (status != algoim_wrapper::CUT)
      {
        const DOMAIN_TYPE facet_domain = (status == algoim_wrapper::UNCUT_NEGATIVE) ? NEG : POS;
        return (facet_domain == dt) ? &SelectIntegrationRule(etfacet, intorder) : nullptr;
      }

      const int npts = wts.size();
      if (npts == 0)
        return nullptr;
      auto ir = new (lh) IntegrationRule(npts, lh);
      for (int i = 0; i < npts; i++)
      {
        Vec<3> x = 0.0;
        for (int d = 0; d < D; d++)
          x(d) = pts[i*D+d];
        double w = wts[i];
        if ((dt == IF) && (D == 2))
        {
          // curve on a 2D facet: physical length element |F J t|, t the unit tangent
          // in facet coordinates, J the Jacobian of the facet-to-element map
          const double gnorm = sqrt(sqr(grads[2*i]) + sqr(grads[2*i+1]));
          if (gnorm == 0)
            w = 0.0;
          else
          {
            Vec<2> t(grads[2*i+1] / gnorm, -grads[2*i] / gnorm);
            IntegrationPoint ip_facet(x, 0.0), ip_vol;
            transform(facetnr, ip_facet, ip_vol);
            FlatMatrix<> J = transform.GetJacobian(facetnr, lh);
            Vec<3> tvol = J * t;
            MappedIntegrationPoint<3,3> mip(ip_vol, trafo);
            w *= L2Norm(mip.GetJacobian() * tvol);
          }
        }
        // (points on 1D facets: weight 1)
        (*ir)[i] = IntegrationPoint(x, w);
      }
      return ir;
    }
  }

  const IntegrationRule *
  AlgoimCutFacetIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                                 const ElementTransformation & trafo,
                                 Facet2ElementTrafo & transform, int facetnr,
                                 DOMAIN_TYPE dt, int intorder, int cf_degree, int subdivlvl,
                                 LocalHeap & lh)
  {
    static Timer timer("AlgoimCutFacetIntegrationRule");
    RegionTimer reg (timer);

    const ELEMENT_TYPE etfacet = transform.FacetType(facetnr);
    if (etfacet == ET_POINT)
    {
      // facets of 1D elements: the point is in the domain or not
      if (dt == IF)
        return nullptr;
      IntegrationPoint ip_facet(0.0, 0.0, 0.0, 0.0), ip_vol;
      transform(facetnr, ip_facet, ip_vol);
      const double val = cf_lset->Evaluate(trafo(ip_vol, lh));
      if ((dt == NEG) != (val < 0))
        return nullptr;
      auto ir = new (lh) IntegrationRule(1, lh);
      (*ir)[0] = IntegrationPoint(0.0, 0.0, 0.0, 1.0);
      return ir;
    }
    const int degree = AlgoimLevelsetDegree(cf_lset, trafo, cf_degree, lh);
    switch (etfacet)
    {
    case ET_SEGM: return T_AlgoimCutFacetIntegrationRule<1>(cf_lset, trafo, transform, facetnr, dt, intorder, degree, subdivlvl, lh);
    case ET_TRIG: case ET_QUAD:
      if (trafo.SpaceDim() != 3)
        throw Exception("AlgoimCutFacetIntegrationRule: 2D facets only for 3D volume elements");
      return T_AlgoimCutFacetIntegrationRule<2>(cf_lset, trafo, transform, facetnr, dt, intorder, degree, subdivlvl, lh);
    default:
      throw Exception("AlgoimCutFacetIntegrationRule: facet type not supported");
    }
  }

  bool AlgoimCutFacetRuleAvailable (ELEMENT_TYPE etfacet)
  {
    return (etfacet == ET_POINT) || (etfacet == ET_SEGM) || (etfacet == ET_TRIG) || (etfacet == ET_QUAD);
  }

  bool AlgoimCutRuleAvailable (const ElementTransformation & trafo, DOMAIN_TYPE dt)
  {
    const ELEMENT_TYPE et = trafo.GetElementType();
    if (!((et == ET_SEGM) || (et == ET_TRIG) || (et == ET_QUAD) || (et == ET_TET) || (et == ET_HEX)))
      return false;
    // interfaces: on volume elements and on boundary elements (co-dimension 2)
    if (dt == IF)
      return trafo.SpaceDim() <= ElementTopology::GetSpaceDim(et) + 1;
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
