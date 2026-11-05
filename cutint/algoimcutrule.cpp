#include "algoimcutrule.hpp"
#include <atomic>
#include "algoim_wrapper.hpp"
#include "../spacetime/SpaceTimeFE.hpp"
#include "../spacetime/SpaceTimeFESpace.hpp"

namespace xintegration
{
  static_assert(NEG == 0 && POS == 1 && IF == 2,
                "algoim_wrapper::CutQuadrature expects NEG = 0, POS = 1, IF = 2");

  void ReportAlgoimUse (bool automatic)
  {
    static std::atomic<bool> reported(false);
    if (reported.exchange(true))
      return;
    cout << IM(0) << "ngsxfem: the quadrature on cut elements uses Algoim by R. I. Saye "
         << "(https://algoim.github.io/)";
    if (automatic)
      cout << IM(0) << ", chosen automatically (cut_quadrature 'auto': for level sets that are not "
           << "piecewise linear GridFunctions and for cases that the standard rules do not support; "
           << "set ngsxfemglobals.cut_quadrature or dCut(..., cut_quadrature=...) to choose explicitly)";
    cout << IM(0) << ". If you publish results, please also cite R. I. Saye, J. Comput. Phys. 448, "
         << "110720, 2022, https://doi.org/10.1016/j.jcp.2021.110720" << endl;
  }

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

    // Number of coordinates that form the simplex of the reference element (triangle, tetrahedron,
    // prism = triangle x interval), see algoim_wrapper::MultiCutQuadrature; 0 for segments,
    // quadrilaterals and hexahedra (the reference element is the unit hypercube).
    inline int SimplexDim (ELEMENT_TYPE et)
    {
      switch (et)
      {
      case ET_TRIG: case ET_PRISM: return 2;
      case ET_TET: return 3;
      default: return 0;
      }
    }

    // polynomial degree in time of the space-time level set function (see AlgoimLevelsetDegree)
    int AlgoimLevelsetTimeDegree (const shared_ptr<CoefficientFunction> & cf_lset, int cf_degree)
    {
      auto gf = dynamic_pointer_cast<GridFunction>(cf_lset);
      if (gf)
      {
        if (auto st_fes = dynamic_cast<SpaceTimeFESpace*>(gf->GetFESpace().get()))
          return max(1, st_fes->GetTimeFE()->Order());
        return 1; // constant in time
      }
      return max(1, cf_degree);
    }

    // Values of the level set function in Algoim's interpolation points of [0,1]^N (for degree):
    // the first DS coordinates are reference coordinates of the element (of its facet facetnr, if
    // transform != nullptr), the last coordinate is the reference time if N > DS. The points are
    // space-time integration points (the time in the weight), at the fixed time tref if N == DS.
    template <int N>
    FlatVector<> LevelsetValuesAtNodes (const shared_ptr<CoefficientFunction> & cf_lset,
                                        const ElementTransformation & trafo,
                                        Facet2ElementTrafo * transform, int facetnr, int DS,
                                        int degree, optional<double> tref, LocalHeap & lh)
    {
      std::vector<double> xi;
      algoim_wrapper::InterpolationPoints(N, degree, xi);
      const int nnodes = xi.size() / N;
      IntegrationRule ir_nodes(nnodes, lh);
      for (int i = 0; i < nnodes; i++)
      {
        Vec<3> p = 0.0;
        for (int d = 0; d < DS; d++)
          p(d) = xi[i*N+d];
        ir_nodes[i] = IntegrationPoint(p, 0.0);
      }
      IntegrationRule & ir_vol = transform ? (*transform)(facetnr, ir_nodes, lh) : ir_nodes;
      // (the facet transformation maps only the first point of a rule on a point facet)
      if (transform && (transform->FacetType(facetnr) == ET_POINT))
        for (int i = 1; i < nnodes; i++)
          ir_vol[i].Point() = ir_vol[0].Point();
      if ((N > DS) || tref)
        for (int i = 0; i < nnodes; i++)
        {
          ir_vol[i].SetWeight((N > DS) ? xi[i*N+N-1] : *tref);
          MarkAsSpaceTimeIntegrationPoint(ir_vol[i]);
        }
      auto & mir_nodes = trafo(ir_vol, lh);
      FlatMatrix<> lsetvals(nnodes, 1, lh);
      cf_lset->Evaluate(mir_nodes, lsetvals);
      return FlatVector<>(nnodes, &lsetvals(0,0));
    }

    // Algoim uses nested one-dimensional Gauss rules (q points) w.r.t. height functions:
    // * for a polynomial of degree intorder on a domain (of dimension dim) bounded by planes, the
    //   integrand of the outermost integral has degree intorder + dim - 1, hence
    //   2q-1 >= intorder + dim - 1,
    // * the (curved) geometry has to be resolved also for low order integrands (e.g. CutInfo
    //   computes cut ratios with intorder 0), hence q > degree.
    inline int AlgoimGaussPoints (int degree, int intorder, int dim)
    {
      return max(degree + 1, (intorder + dim + 1) / 2);
    }

    // Space-time rule on the element (spatial dimension D) times the reference time interval
    // [0,1], with the conventions of SpaceTimeCutIntegrationRule: spatial points, the reference
    // time in the weight of the integration point (marked as space-time point), weights in the
    // array. Interfaces: the measure is the iterated one (int_0^1 int_{Gamma(t)} ds dt), the
    // space-time surface weights w.r.t. (x,t) are multiplied by |grad_x phi| / |grad_(x,t) phi|
    // (coarea formula) and the spatial transformation factor |F^{-T} n_x|.
    template <int D>
    tuple<const IntegrationRule *, Array<double>>
    T_AlgoimSpaceTimeCutIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                                         const ElementTransformation & trafo,
                                         DOMAIN_TYPE dt, int intorder, int time_intorder,
                                         int degree, int subdivlvl, LocalHeap & lh)
    {
      constexpr int N = D + 1;
      const ELEMENT_TYPE et = trafo.GetElementType();
      const int sdim = SimplexDim(et);
      FlatVector<> lsetvals = LevelsetValuesAtNodes<N>(cf_lset, trafo, nullptr, 0, D, degree, nullopt, lh);

      const int q = AlgoimGaussPoints(degree, max(intorder, time_intorder), N);
      std::vector<double> pts, wts, grads;
      const double * valptr = &lsetvals(0);
      const int domain = int(dt);
      auto status = algoim_wrapper::MultiCutQuadrature(N, 1, &degree, &valptr, &domain, sdim, q, subdivlvl,
                                                       pts, wts, grads, true);
      if (status == algoim_wrapper::EMPTY)
        return make_tuple(nullptr, Array<double>());

      if (status == algoim_wrapper::FULL)
      {
        // tensor product of the standard rules in space and time
        const IntegrationRule & ir_space = SelectIntegrationRule(et, intorder);
        const IntegrationRule & ir_time = SelectIntegrationRule(ET_SEGM, time_intorder);
        auto ir = new (lh) IntegrationRule(ir_space.Size() * ir_time.Size(), lh);
        Array<double> wei_arr(ir->Size());
        int cnt = 0;
        for (auto & ipt : ir_time)
          for (auto & ips : ir_space)
          {
            (*ir)[cnt] = IntegrationPoint(ips.Point(), ipt(0));
            MarkAsSpaceTimeIntegrationPoint((*ir)[cnt]);
            wei_arr[cnt] = ips.Weight() * ipt.Weight();
            cnt++;
          }
        return make_tuple(ir, std::move(wei_arr));
      }

      const int npts = wts.size();
      if (npts == 0)
        return make_tuple(nullptr, Array<double>());
      auto ir = new (lh) IntegrationRule(npts, lh);
      Array<double> wei_arr(npts);
      for (int i = 0; i < npts; i++)
      {
        Vec<3> x = 0.0;
        for (int d = 0; d < D; d++)
          x(d) = pts[i*N+d];
        const double t = pts[i*N+D];
        double w = wts[i];
        IntegrationPoint ip(x, t);
        MarkAsSpaceTimeIntegrationPoint(ip);
        if (dt == IF)
        {
          Vec<D> gx;
          double g2 = sqr(grads[i*N+D]);
          for (int d = 0; d < D; d++)
          {
            gx(d) = grads[i*N+d];
            g2 += sqr(gx(d));
          }
          // |grad_x phi| / |grad phi| * |F^{-T} grad_x phi| / |grad_x phi|
          const double gxnorm = L2Norm(gx);
          w = (gxnorm > 0) ? w * gxnorm / sqrt(g2) * InterfaceWeightFactor<D,D>(ip, trafo, gx) : 0.0;
        }
        (*ir)[i] = ip;
        wei_arr[i] = w;
      }
      return make_tuple(ir, std::move(wei_arr));
    }

    template <int D>
    tuple<const IntegrationRule *, Array<double>>
    T_AlgoimCutIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                                const ElementTransformation & trafo,
                                DOMAIN_TYPE dt,
                                int intorder,
                                int degree,
                                int subdivlvl,
                                optional<double> tref,
                                LocalHeap & lh)
    {
      const ELEMENT_TYPE et = trafo.GetElementType();
      const int sdim = SimplexDim(et);

      // level set values in the interpolation points in [0,1]^D (for simplices,
      // these points are also outside of the reference element, where the
      // level set function is evaluated by polynomial extension); at the
      // reference time tref (space-time level sets), if given
      FlatVector<> lsetvals = LevelsetValuesAtNodes<D>(cf_lset, trafo, nullptr, 0, D, degree, tref, lh);

      // Algoim's quadrature on the reference element
      const int q = AlgoimGaussPoints(degree, intorder, D);
      std::vector<double> pts, wts, grads;
      const double * valptr = &lsetvals(0);
      const int domain = int(dt);
      auto status = algoim_wrapper::MultiCutQuadrature(D, 1, &degree, &valptr, &domain, sdim, q, subdivlvl,
                                                       pts, wts, grads);

      if (status != algoim_wrapper::PARTIAL)
      {
        if (status == algoim_wrapper::EMPTY)
          return make_tuple(nullptr, Array<double>());
        const IntegrationRule & ir = SelectIntegrationRule(et, intorder);
        Array<double> wei_arr(ir.Size());
        for (int i = 0; i < ir.Size(); i++)
          wei_arr[i] = ir[i].Weight();
        if (!tref)
          return make_tuple(&ir, std::move(wei_arr));
        // fixed reference time: copy of the rule with the time in the weights
        auto ir_t = new (lh) IntegrationRule(ir.Size(), lh);
        for (int i = 0; i < ir.Size(); i++)
        {
          (*ir_t)[i] = IntegrationPoint(ir[i].Point(), *tref);
          MarkAsSpaceTimeIntegrationPoint((*ir_t)[i]);
        }
        return make_tuple(ir_t, std::move(wei_arr));
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
          IntegrationPoint ip(x, tref ? *tref : 0.0);
          if (tref)
            MarkAsSpaceTimeIntegrationPoint(ip);
          if (trafo.SpaceDim() == D)
            w *= InterfaceWeightFactor<D,D>(ip, trafo, g);
          else if constexpr (D < 3)
            w *= InterfaceWeightFactor<D,D+1>(ip, trafo, g);
        }
        // with a fixed reference time, the weight of the point is the time (weights in the array)
        (*ir)[i] = IntegrationPoint(x, tref ? *tref : w);
        if (tref)
          MarkAsSpaceTimeIntegrationPoint((*ir)[i]);
        wei_arr[i] = w;
      }
      return make_tuple(ir, std::move(wei_arr));
    }

    // Time slicing (for space-time cells of dimension 4: Algoim's quadrature in dimension 4 w.r.t.
    // several polynomials, e.g. a level set and the simplex constraint of tetrahedra, is very
    // expensive and not robust, the elimination steps lead to resultants of high degree with many
    // roots): the time interval [0,1] is split at the times where one of the level set functions
    // changes its sign in a vertex of the element (where the topology of the cut changes, as for
    // the standard space-time rules); on each sub-interval, the Gauss rule in time (order
    // time_intorder) is combined with the spatial rules at the fixed times (spatial_rule(t, lh),
    // with the conventions of the space-time rules). For interfaces, the spatial rules directly
    // give the iterated measure int_0^1 int_{Gamma(t)} ds dt.
    Array<double> TimeBreakpoints (FlatArray<shared_ptr<CoefficientFunction>> lsets, FlatArray<int> degrees_time,
                                   const ElementTransformation & trafo, LocalHeap & lh)
    {
      HeapReset hr(lh);
      const ELEMENT_TYPE et = trafo.GetElementType();
      const int nv = ElementTopology::GetNVertices(et);
      const POINT3D * verts = ElementTopology::GetVertices(et);
      Array<double> tbreak;
      tbreak.Append(0.0);
      for (int l = 0; l < lsets.Size(); l++)
      {
        std::vector<double> tnodes;
        algoim_wrapper::InterpolationPoints(1, degrees_time[l], tnodes);
        const int nt = tnodes.size();
        IntegrationRule ir_vt(nv * nt, lh);
        for (int v = 0; v < nv; v++)
          for (int j = 0; j < nt; j++)
          {
            ir_vt[v*nt+j] = IntegrationPoint(verts[v][0], verts[v][1], verts[v][2], tnodes[j]);
            MarkAsSpaceTimeIntegrationPoint(ir_vt[v*nt+j]);
          }
        auto & mir_vt = trafo(ir_vt, lh);
        FlatMatrix<> vals(nv * nt, 1, lh);
        lsets[l]->Evaluate(mir_vt, vals);
        for (int v = 0; v < nv; v++)
        {
          std::vector<double> roots, wts, grads;
          auto status = algoim_wrapper::CutQuadrature(1, degrees_time[l], &vals(v*nt,0), false, int(IF), 1, 0,
                                                      roots, wts, grads);
          if (status == algoim_wrapper::CUT)
            for (double r : roots)
              tbreak.Append(r);
        }
      }
      tbreak.Append(1.0);
      QuickSort(tbreak);
      return tbreak;
    }

    template <typename SPATIAL_RULE>
    tuple<const IntegrationRule *, Array<double>>
    TimeSlicedRule (FlatArray<double> tbreak, int time_intorder, SPATIAL_RULE && spatial_rule, LocalHeap & lh)
    {
      const IntegrationRule & ir_time = SelectIntegrationRule(ET_SEGM, time_intorder);
      Array<IntegrationPoint> points;
      Array<double> wei_arr;
      for (int k = 0; k + 1 < tbreak.Size(); k++)
      {
        const double t0 = tbreak[k], t1 = tbreak[k+1];
        if (t1 - t0 < 1e-14)
          continue;
        for (auto & ipt : ir_time)
        {
          HeapReset hr(lh);
          const double t = t0 + (t1 - t0) * ipt(0);
          const double wt = (t1 - t0) * ipt.Weight();
          auto [ir_t, wei_t] = spatial_rule(t, lh);
          if (ir_t == nullptr)
            continue;
          for (int i = 0; i < ir_t->Size(); i++)
          {
            points.Append((*ir_t)[i]); // spatial point, time t in the weight (marked)
            wei_arr.Append(wei_t[i] * wt);
          }
        }
      }
      if (points.Size() == 0)
        return make_tuple(nullptr, Array<double>());
      auto ir = new (lh) IntegrationRule(points.Size(), lh);
      for (int i = 0; i < points.Size(); i++)
        (*ir)[i] = points[i];
      return make_tuple(ir, std::move(wei_arr));
    }

    // space-time rule of a single level set by time slicing (tetrahedra, prisms)
    template <int D>
    tuple<const IntegrationRule *, Array<double>>
    T_AlgoimSpaceTimeCutIntegrationRuleTimeSlices (const shared_ptr<CoefficientFunction> & cf_lset,
                                                   const ElementTransformation & trafo,
                                                   DOMAIN_TYPE dt, int intorder, int time_intorder,
                                                   int degree_space, int degree_time, int subdivlvl,
                                                   LocalHeap & lh)
    {
      Array<double> tbreak = TimeBreakpoints(FlatArray<shared_ptr<CoefficientFunction>>(1, const_cast<shared_ptr<CoefficientFunction>*>(&cf_lset)),
                                             FlatArray<int>(1, &degree_time), trafo, lh);
      return TimeSlicedRule(tbreak, time_intorder,
                            [&] (double t, LocalHeap & lh)
                            {
                              return T_AlgoimCutIntegrationRule<D>(cf_lset, trafo, dt, intorder, degree_space,
                                                                   subdivlvl, t, lh);
                            }, lh);
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
      FlatVector<> lsetvals = LevelsetValuesAtNodes<D>(cf_lset, trafo, &transform, facetnr, D, degree, nullopt, lh);

      const int q = AlgoimGaussPoints(degree, intorder, D);
      std::vector<double> pts, wts, grads;
      auto status = algoim_wrapper::CutQuadrature(D, degree, &lsetvals(0), simplex, int(dt), q,
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

    // Space-time cut rule on the facet facetnr (of dimension DF) of the element times the
    // reference time interval, with the conventions of SpaceTimeCutIntegrationRuleUntransformed:
    // points in reference coordinates of the facet with the reference time in the weight of the
    // integration point, weights in the array:
    //  * NEG/POS: weights w.r.t. the reference facet (times the time weight),
    //  * IF: weights w.r.t. the iterated physical measure int_0^1 int_{Gamma(t) on the facet} dt,
    //    i.e. the weight of the space-time curve / surface times |grad_x phi| / |grad_(x,t) phi|
    //    times the physical measure of the spatial interface (points on 1D facets: 1, curves
    //    on 2D facets: |F J t|, t the unit tangent in facet coordinates).
    template <int DF>
    tuple<const IntegrationRule *, Array<double>>
    T_AlgoimSpaceTimeCutFacetIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                                              const ElementTransformation & trafo,
                                              Facet2ElementTrafo & transform, int facetnr,
                                              DOMAIN_TYPE dt, int intorder, int time_intorder,
                                              int degree, int subdivlvl, LocalHeap & lh)
    {
      constexpr int N = DF + 1;
      const ELEMENT_TYPE etfacet = transform.FacetType(facetnr);
      const int sdim = (etfacet == ET_TRIG) ? 2 : 0;
      FlatVector<> lsetvals = LevelsetValuesAtNodes<N>(cf_lset, trafo, &transform, facetnr, DF, degree,
                                                       nullopt, lh);

      const int q = AlgoimGaussPoints(degree, max(intorder, time_intorder), N);
      std::vector<double> pts, wts, grads;
      const double * valptr = &lsetvals(0);
      const int domain = int(dt);
      auto status = algoim_wrapper::MultiCutQuadrature(N, 1, &degree, &valptr, &domain, sdim, q, subdivlvl,
                                                       pts, wts, grads, true);
      if (status == algoim_wrapper::EMPTY)
        return make_tuple(nullptr, Array<double>());

      if (status == algoim_wrapper::FULL)
      {
        // tensor product of the standard rules on the facet and in time
        IntegrationRule ir_point(1, lh);
        ir_point[0] = IntegrationPoint(0.0, 0.0, 0.0, 1.0);
        const IntegrationRule & ir_space = (DF == 0) ? ir_point : SelectIntegrationRule(etfacet, intorder);
        const IntegrationRule & ir_time = SelectIntegrationRule(ET_SEGM, time_intorder);
        auto ir = new (lh) IntegrationRule(ir_space.Size() * ir_time.Size(), lh);
        Array<double> wei_arr(ir->Size());
        int cnt = 0;
        for (auto & ipt : ir_time)
          for (auto & ips : ir_space)
          {
            (*ir)[cnt] = IntegrationPoint(ips.Point(), ipt(0));
            MarkAsSpaceTimeIntegrationPoint((*ir)[cnt]);
            wei_arr[cnt] = ips.Weight() * ipt.Weight();
            cnt++;
          }
        return make_tuple(ir, std::move(wei_arr));
      }

      const int npts = wts.size();
      if (npts == 0)
        return make_tuple(nullptr, Array<double>());
      auto ir = new (lh) IntegrationRule(npts, lh);
      Array<double> wei_arr(npts);
      for (int i = 0; i < npts; i++)
      {
        Vec<3> x = 0.0;
        for (int d = 0; d < DF; d++)
          x(d) = pts[i*N+d];
        const double t = pts[i*N+DF];
        double w = wts[i];
        if ((dt == IF) && (DF > 0))
        {
          Vec<DF> gx;
          double g2 = sqr(grads[i*N+DF]);
          for (int d = 0; d < DF; d++)
          {
            gx(d) = grads[i*N+d];
            g2 += sqr(gx(d));
          }
          const double gxnorm = L2Norm(gx);
          w = (gxnorm > 0) ? w * gxnorm / sqrt(g2) : 0.0;
          if constexpr (DF == 2)
            if (gxnorm > 0)
            {
              // curve on a 2D facet: physical length element |F J t|
              Vec<2> tang(gx(1) / gxnorm, -gx(0) / gxnorm);
              IntegrationPoint ip_facet(x, t), ip_vol;
              transform(facetnr, ip_facet, ip_vol);
              ip_vol.SetWeight(t);
              MarkAsSpaceTimeIntegrationPoint(ip_vol);
              FlatMatrix<> J = transform.GetJacobian(facetnr, lh);
              Vec<3> tvol = J * tang;
              MappedIntegrationPoint<3,3> mip(ip_vol, trafo);
              w *= L2Norm(mip.GetJacobian() * tvol);
            }
        }
        (*ir)[i] = IntegrationPoint(x, t);
        MarkAsSpaceTimeIntegrationPoint((*ir)[i]);
        wei_arr[i] = w;
      }
      return make_tuple(ir, std::move(wei_arr));
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
    ReportAlgoimUse();

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

  tuple<const IntegrationRule *, Array<double>>
  AlgoimSpaceTimeCutFacetIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                                          const ElementTransformation & trafo,
                                          Facet2ElementTrafo & transform, int facetnr,
                                          DOMAIN_TYPE dt, int intorder, int time_intorder,
                                          int cf_degree, int subdivlvl, LocalHeap & lh)
  {
    static Timer timer("AlgoimSpaceTimeCutFacetIntegrationRule");
    RegionTimer reg (timer);
    ReportAlgoimUse();

    const ELEMENT_TYPE etfacet = transform.FacetType(facetnr);
    // the interface meets a point facet only at isolated times (zero iterated measure)
    if ((etfacet == ET_POINT) && (dt == IF))
      return make_tuple(nullptr, Array<double>());
    const int degree = AlgoimLevelsetDegree(cf_lset, trafo, cf_degree, lh, true);
    switch (etfacet)
    {
    case ET_POINT: return T_AlgoimSpaceTimeCutFacetIntegrationRule<0>(cf_lset, trafo, transform, facetnr, dt, intorder, time_intorder, degree, subdivlvl, lh);
    case ET_SEGM: return T_AlgoimSpaceTimeCutFacetIntegrationRule<1>(cf_lset, trafo, transform, facetnr, dt, intorder, time_intorder, degree, subdivlvl, lh);
    case ET_TRIG: case ET_QUAD:
      if (trafo.SpaceDim() != 3)
        throw Exception("AlgoimSpaceTimeCutFacetIntegrationRule: 2D facets only for 3D volume elements");
      return T_AlgoimSpaceTimeCutFacetIntegrationRule<2>(cf_lset, trafo, transform, facetnr, dt, intorder, time_intorder, degree, subdivlvl, lh);
    default:
      throw Exception("AlgoimSpaceTimeCutFacetIntegrationRule: facet type not supported");
    }
  }

  namespace
  {
    // Rule for the union of the parts given by the domain type tuples of lsetintdom (multiple
    // level sets): spatial (N == D; at the fixed reference time tref, if given) or space-time
    // (N == D + 1, the time in the last coordinate of Algoim's cell), with the conventions of
    // AlgoimCutIntegrationRule / AlgoimSpaceTimeCutIntegrationRule.
    template <int D, int N>
    tuple<const IntegrationRule *, Array<double>>
    T_AlgoimMultiLevelsetRule (const LevelsetIntegrationDomain & lsetintdom,
                               const ElementTransformation & trafo,
                               int cf_degree, optional<double> tref, LocalHeap & lh)
    {
      static_assert((N == D) || (N == D + 1), "spatial or space-time rule");
      constexpr bool spacetime = (N > D);
      const bool timepoints = spacetime || tref.has_value(); // the time in the weights of the points
      const ELEMENT_TYPE et = trafo.GetElementType();
      const int sdim = SimplexDim(et);
      const auto & gflsets = lsetintdom.GetLevelsetGFs();
      const int M = gflsets.Size();
      const int intorder = max(lsetintdom.GetIntegrationOrder(), 0);
      const int time_intorder = max(lsetintdom.GetTimeIntegrationOrder(), 0);

      // level set values in the interpolation points (per level set, depending on its degree)
      Array<int> degrees(M);
      std::vector<std::vector<double>> values(M);
      Array<const double*> valptrs(M);
      int maxdegree = 1;
      for (int l = 0; l < M; l++)
      {
        degrees[l] = AlgoimLevelsetDegree(gflsets[l], trafo, cf_degree, lh, spacetime);
        maxdegree = max(maxdegree, degrees[l]);
        FlatVector<> lsetvals = LevelsetValuesAtNodes<N>(gflsets[l], trafo, nullptr, 0, D, degrees[l],
                                                         spacetime ? nullopt : tref, lh);
        values[l].assign(&lsetvals(0), &lsetvals(0) + lsetvals.Size());
        valptrs[l] = values[l].data();
      }

      // rules for all parts (domain type tuples), the union is returned
      const int q = AlgoimGaussPoints(maxdegree, spacetime ? max(intorder, time_intorder) : intorder, N);
      Array<IntegrationPoint> points;
      Array<double> weights;
      auto append = [&] (const IntegrationPoint & ip, double w)
      {
        points.Append(ip);
        if (timepoints)
          MarkAsSpaceTimeIntegrationPoint(points.Last());
        weights.Append(w);
      };
      for (const Array<DOMAIN_TYPE> & dts : lsetintdom.GetDomainTypes())
      {
        Array<int> domains(M);
        int iface = -1, ninterfaces = 0;
        for (int l = 0; l < M; l++)
        {
          domains[l] = int(dts[l]);
          if (dts[l] == IF)
          {
            iface = l;
            ninterfaces++;
          }
        }
        if (ninterfaces > 1)
          throw Exception("Algoim-based cut rules for multiple level sets: domains with more than one "
                          "interface condition (co-dimension > 1) are not supported (only for piecewise "
                          "linear level set GridFunctions with cut_quadrature = 'standard')");

        std::vector<double> pts, wts, grads;
        auto status = algoim_wrapper::MultiCutQuadrature(N, M, degrees.Data(), valptrs.Data(), domains.Data(),
                                                         sdim, q, lsetintdom.GetNSubdivisionLevels(),
                                                         pts, wts, grads, spacetime);
        if (status == algoim_wrapper::EMPTY)
          continue;
        if (status == algoim_wrapper::FULL)
        {
          const IntegrationRule & ir_space = SelectIntegrationRule(et, intorder);
          if constexpr (spacetime)
          {
            for (auto & ipt : SelectIntegrationRule(ET_SEGM, time_intorder))
              for (auto & ips : ir_space)
                append(IntegrationPoint(ips.Point(), ipt(0)), ips.Weight() * ipt.Weight());
          }
          else
            for (auto & ips : ir_space)
              append(IntegrationPoint(ips.Point(), tref ? *tref : ips.Weight()), ips.Weight());
          continue;
        }
        for (size_t i = 0; i < wts.size(); i++)
        {
          Vec<3> x = 0.0;
          for (int d = 0; d < D; d++)
            x(d) = pts[i*N+d];
          const double t = spacetime ? pts[i*N+D] : (tref ? *tref : 0.0);
          double w = wts[i];
          if (iface >= 0)
          {
            Vec<D> gx;
            for (int d = 0; d < D; d++)
              gx(d) = grads[i*N+d];
            IntegrationPoint ip(x, t);
            if (timepoints)
              MarkAsSpaceTimeIntegrationPoint(ip);
            if constexpr (spacetime)
            {
              // coarea factor |grad_x phi| / |grad_(x,t) phi| (cf. T_AlgoimSpaceTimeCutIntegrationRule)
              const double gxnorm = L2Norm(gx);
              const double gnorm = sqrt(sqr(gxnorm) + sqr(grads[i*N+D]));
              w = (gxnorm > 0) ? w * gxnorm / gnorm * InterfaceWeightFactor<D,D>(ip, trafo, gx) : 0.0;
            }
            else if (trafo.SpaceDim() == D)
              w *= InterfaceWeightFactor<D,D>(ip, trafo, gx);
            else if constexpr (D < 3)
              w *= InterfaceWeightFactor<D,D+1>(ip, trafo, gx);
          }
          append(IntegrationPoint(x, timepoints ? t : w), w);
        }
      }

      if (points.Size() == 0)
        return make_tuple(nullptr, Array<double>());
      auto ir = new (lh) IntegrationRule(points.Size(), lh);
      for (int i = 0; i < points.Size(); i++)
        (*ir)[i] = points[i];
      return make_tuple(ir, std::move(weights));
    }
  }

  tuple<const IntegrationRule *, Array<double>>
  AlgoimMultiLevelsetCutIntegrationRule (const LevelsetIntegrationDomain & lsetintdom,
                                         const ElementTransformation & trafo,
                                         int cf_degree, LocalHeap & lh)
  {
    static Timer timer("AlgoimMultiLevelsetCutIntegrationRule");
    RegionTimer reg (timer);
    ReportAlgoimUse();
    const int D = ElementTopology::GetSpaceDim(trafo.GetElementType());
    if (lsetintdom.GetTimeIntegrationOrder() >= 0)
    {
      if (lsetintdom.HasReferenceTime())
        throw Exception("space-time quadrature rule shall not have a fixed reference time.");
      switch (D)
      {
      case 1: return T_AlgoimMultiLevelsetRule<1,2>(lsetintdom, trafo, cf_degree, nullopt, lh);
      case 2: return T_AlgoimMultiLevelsetRule<2,3>(lsetintdom, trafo, cf_degree, nullopt, lh);
      case 3:
      {
        // time slicing (see TimeBreakpoints)
        const auto & gflsets = lsetintdom.GetLevelsetGFs();
        Array<shared_ptr<CoefficientFunction>> lsets(gflsets.Size());
        Array<int> degrees_time(gflsets.Size());
        for (int l = 0; l < gflsets.Size(); l++)
        {
          lsets[l] = gflsets[l];
          degrees_time[l] = AlgoimLevelsetTimeDegree(gflsets[l], cf_degree);
        }
        Array<double> tbreak = TimeBreakpoints(lsets, degrees_time, trafo, lh);
        return TimeSlicedRule(tbreak, lsetintdom.GetTimeIntegrationOrder(),
                              [&] (double t, LocalHeap & lh)
                              {
                                return T_AlgoimMultiLevelsetRule<3,3>(lsetintdom, trafo, cf_degree, t, lh);
                              }, lh);
      }
      default:
        throw Exception("AlgoimMultiLevelsetCutIntegrationRule: dimension not supported");
      }
    }
    optional<double> tref;
    if (lsetintdom.HasReferenceTime())
      tref = lsetintdom.ReferenceTime();
    switch (D)
    {
    case 1: return T_AlgoimMultiLevelsetRule<1,1>(lsetintdom, trafo, cf_degree, tref, lh);
    case 2: return T_AlgoimMultiLevelsetRule<2,2>(lsetintdom, trafo, cf_degree, tref, lh);
    case 3: return T_AlgoimMultiLevelsetRule<3,3>(lsetintdom, trafo, cf_degree, tref, lh);
    default:
      throw Exception("AlgoimMultiLevelsetCutIntegrationRule: dimension not supported");
    }
  }

  bool AlgoimCutFacetRuleAvailable (ELEMENT_TYPE etfacet)
  {
    return (etfacet == ET_POINT) || (etfacet == ET_SEGM) || (etfacet == ET_TRIG) || (etfacet == ET_QUAD);
  }

  bool AlgoimCutRuleAvailable (const ElementTransformation & trafo, DOMAIN_TYPE dt)
  {
    const ELEMENT_TYPE et = trafo.GetElementType();
    if (!((et == ET_SEGM) || (et == ET_TRIG) || (et == ET_QUAD) || (et == ET_TET) || (et == ET_HEX)
          || (et == ET_PRISM)))
      return false;
    // interfaces: on volume elements and on boundary elements (co-dimension 2)
    if (dt == IF)
      return trafo.SpaceDim() <= ElementTopology::GetSpaceDim(et) + 1;
    return true;
  }

  bool AlgoimSpaceTimeCutRuleAvailable (const ElementTransformation & trafo, DOMAIN_TYPE dt)
  {
    if (!AlgoimCutRuleAvailable(trafo, dt))
      return false;
    // space-time interfaces only on volume elements
    return (dt != IF) || (trafo.SpaceDim() == ElementTopology::GetSpaceDim(trafo.GetElementType()));
  }

  int AlgoimLevelsetDegree (const shared_ptr<CoefficientFunction> & cf_lset,
                            const ElementTransformation & trafo,
                            int cf_degree,
                            LocalHeap & lh,
                            bool spacetime)
  {
    auto gf = dynamic_pointer_cast<GridFunction>(cf_lset);
    if (gf && (gf->Dimension() == 1))
    {
      if (auto st_fes = dynamic_cast<SpaceTimeFESpace*>(gf->GetFESpace().get()))
      {
        int order = st_fes->GetSpaceFESpace()->GetFE(trafo.GetElementId(), lh).Order();
        if (spacetime)
          order = max(order, st_fes->GetTimeFE()->Order());
        return max(1, order);
      }
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
                            LocalHeap & lh,
                            optional<double> tref)
  {
    static Timer timer("AlgoimCutIntegrationRule");
    RegionTimer reg (timer);
    ReportAlgoimUse();

    if (!AlgoimCutRuleAvailable(trafo, dt))
      throw Exception("AlgoimCutIntegrationRule: element type / domain type combination not supported");

    const int degree = AlgoimLevelsetDegree(cf_lset, trafo, cf_degree, lh);
    switch (ElementTopology::GetSpaceDim(trafo.GetElementType()))
    {
    case 1: return T_AlgoimCutIntegrationRule<1>(cf_lset, trafo, dt, intorder, degree, subdivlvl, tref, lh);
    case 2: return T_AlgoimCutIntegrationRule<2>(cf_lset, trafo, dt, intorder, degree, subdivlvl, tref, lh);
    case 3: return T_AlgoimCutIntegrationRule<3>(cf_lset, trafo, dt, intorder, degree, subdivlvl, tref, lh);
    default:
      throw Exception("AlgoimCutIntegrationRule: dimension not supported");
    }
  }

  tuple<const IntegrationRule *, Array<double>>
  AlgoimSpaceTimeCutIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                                     const ElementTransformation & trafo,
                                     DOMAIN_TYPE dt,
                                     int intorder,
                                     int time_intorder,
                                     int cf_degree,
                                     int subdivlvl,
                                     LocalHeap & lh)
  {
    static Timer timer("AlgoimSpaceTimeCutIntegrationRule");
    RegionTimer reg (timer);
    ReportAlgoimUse();

    if (!AlgoimSpaceTimeCutRuleAvailable(trafo, dt))
      throw Exception("AlgoimSpaceTimeCutIntegrationRule: element type / domain type combination not supported");

    // tetrahedra and prisms: time slicing (see T_AlgoimSpaceTimeCutIntegrationRuleTimeSlices)
    if ((trafo.GetElementType() == ET_TET) || (trafo.GetElementType() == ET_PRISM))
      return T_AlgoimSpaceTimeCutIntegrationRuleTimeSlices<3>(cf_lset, trafo, dt, intorder, time_intorder,
                                                              AlgoimLevelsetDegree(cf_lset, trafo, cf_degree, lh),
                                                              AlgoimLevelsetTimeDegree(cf_lset, cf_degree),
                                                              subdivlvl, lh);
    const int degree = AlgoimLevelsetDegree(cf_lset, trafo, cf_degree, lh, true);
    switch (ElementTopology::GetSpaceDim(trafo.GetElementType()))
    {
    case 1: return T_AlgoimSpaceTimeCutIntegrationRule<1>(cf_lset, trafo, dt, intorder, time_intorder, degree, subdivlvl, lh);
    case 2: return T_AlgoimSpaceTimeCutIntegrationRule<2>(cf_lset, trafo, dt, intorder, time_intorder, degree, subdivlvl, lh);
    case 3: return T_AlgoimSpaceTimeCutIntegrationRule<3>(cf_lset, trafo, dt, intorder, time_intorder, degree, subdivlvl, lh);
    default:
      throw Exception("AlgoimSpaceTimeCutIntegrationRule: dimension not supported");
    }
  }
}
