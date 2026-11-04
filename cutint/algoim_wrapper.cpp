/*
  The only translation unit of ngsxfem that includes the bundled Algoim library
  (some Algoim headers define non-inline functions). See algoim_wrapper.hpp.

  Algoim is due to R. I. Saye, see https://algoim.github.io/ and
  external_dependencies/algoim/README.md (license, citation).
*/

#include "algoim_wrapper.hpp"

// Algoim includes <lapacke.h>, which has to resolve to cutint/algoim_lapacke/lapacke.h
#include "../external_dependencies/algoim/algoim/quadrature_multipoly.hpp"

#ifndef NGSXFEM_ALGOIM_LAPACKE
#error "Algoim picked up a lapacke.h other than cutint/algoim_lapacke/lapacke.h; check the include directories of xfem_core."
#endif

#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>

namespace xintegration
{
  namespace algoim_wrapper
  {
    namespace
    {
      using algoim::real;
      using algoim::uvector;
      using algoim::xarray;
      using algoim::MultiLoop;
      using algoim::SparkStack; // needed by the macro algoim_spark_alloc

      template <int N>
      void T_InterpolationPoints (int degree, std::vector<double> & points)
      {
        const int P = degree + 1;
        int npoints = 1;
        for (int d = 0; d < N; d++)
          npoints *= P;
        points.resize(npoints * N);
        int k = 0;
        for (MultiLoop<N> i(0, P); ~i; ++i, ++k)
          for (int d = 0; d < N; d++)
            points[k*N+d] = algoim::bernstein::modifiedChebyshevNode(i(d), P);
      }

      // relative size of the generic perturbation of the level sets for several level sets,
      // see PerturbGenerically (the geometry is perturbed by this relative amount)
      const real multi_perturbation = 1e-9;

      // Gauss-Legendre quadrature on all one-dimensional sub-intervals: Algoim's alternative
      // "AutoMixed" (tanh-sinh where vertical tangents are detected) pays off for large q only;
      // for the moderate q of finite element computations it was less accurate than "AlwaysGL"
      // in our tests (2D and 3D, all element types)
      const algoim::QuadStrategy strategy = algoim::AlwaysGL;

      // true if Algoim detected (possible) vertical tangents of the interface in the box, i.e. a
      // discriminant with non-empty mask, w.r.t. any coordinate direction (on the top level) or
      // w.r.t. a height direction chosen on a lower level. If vertical tangents are present for
      // a chosen height direction, the integrands of the base integrals have square-root type
      // singularities and the Gauss rules converge only algebraically w.r.t. q. If they are
      // present for other directions, the interface is strongly curved relative to the box,
      // which typically comes with nearly vertical tangents for the chosen direction.
      template <int N>
      bool HasVerticalTangentsInHierarchy (const algoim::ImplicitPolyQuadrature<N> & ipquad)
      {
        if constexpr (N == 1)
          return false;
        else
        {
          if (ipquad.k < 0 || ipquad.k == N)
            return false;
          return ipquad.base.auto_apply_TS || HasVerticalTangentsInHierarchy<N-1>(ipquad.base);
        }
      }

      template <int N>
      bool HasVerticalTangents (algoim::ImplicitPolyQuadrature<N> & ipquad)
      {
        if constexpr (N == 1)
          return false;
        else
        {
          if (ipquad.k < 0 || ipquad.k == N)
            return false;
          uvector<bool,N> has_disc;
          algoim::detail::score_estimate(ipquad.phi, has_disc);
          return any(has_disc) || HasVerticalTangentsInHierarchy<N>(ipquad);
        }
      }

      // Reduces the degree of phi (per direction, up to a small tolerance) as long as it stays at
      // least 1 (as algoim::bernstein::autoReduction, which may also reduce to degree 0; Algoim's
      // elimination steps do not handle polynomials that are constant in a direction correctly,
      // e.g. the level set y - c of a horizontal line).
      template <int N>
      void ReduceDegree (xarray<real,N> & alpha, real tol = 1.0e3 * std::numeric_limits<real>::epsilon())
      {
        for (int dim = 0; dim < N; dim++)
          while (alpha.ext(dim) > 2)
          {
            xarray<real,N> beta(nullptr, alpha.ext()), gamma(nullptr, alpha.ext());
            algoim_spark_alloc(real, beta, gamma);
            beta = alpha;
            algoim::bernstein::bernsteinReduction(beta, dim);
            algoim::bernstein::bernsteinElevate(beta, gamma);
            gamma -= alpha;
            if (std::sqrt(std::abs(algoim::bernstein::squaredL2norm(gamma))) >= tol * std::sqrt(std::abs(algoim::bernstein::squaredL2norm(alpha))))
              break;
            alpha.alterExtent(beta.ext());
            alpha = beta;
          }
      }

      // As ImplicitPolyQuadrature::build(true, false) (outer integral), but the height direction
      // is chosen w.r.t. the interface polynomial p only (with Algoim's heuristic), not w.r.t. all
      // polynomials of the set. The other polynomials only restrict the interface {p = 0}; Algoim's
      // combined choice can result in a height direction in which p is constant (e.g. a
      // horizontal line), which is not detected as degenerate and loses parts of the interface.
      template <int N>
      void BuildForInterface (algoim::ImplicitPolyQuadrature<N> & ipquad,
                              const xarray<real,N> & p, const algoim::booluarray<N,ALGOIM_M> & mask)
      {
        using IPQ = algoim::ImplicitPolyQuadrature<N>;
        if constexpr (N == 1)
          ipquad.build(true, false);
        else
        {
          algoim::PolySet<N,ALGOIM_M> pset;
          pset.push_back(p, mask);
          uvector<bool,N> has_disc_p, has_disc;
          uvector<real,N> score = algoim::detail::score_estimate(pset, has_disc_p);
          score /= 2 * max(abs(score));
          for (int i = 0; i < N; ++i)
            if (!has_disc_p(i))
              score(i) += 1.0;
          const int k = argmax(score);

          // the remainder of build(true, false) with this height direction; the dimension-
          // aggregated surface scheme is only needed if {p = 0} has vertical tangents (vertical
          // tangents of the other polynomials are resolved by the base integrals)
          algoim::detail::score_estimate(ipquad.phi, has_disc);
          ipquad.type = IPQ::OuterSingle;
          ipquad.auto_apply_TS = false;
          ipquad.k = k;
          algoim::detail::eliminate_axis(ipquad.phi, k, ipquad.base.phi);
          ipquad.base.build(false, has_disc(k));
          if (has_disc_p(k))
          {
            ipquad.type = IPQ::OuterAggregate;
            for (int i = 0; i < N; ++i)
              if (i != k)
              {
                auto & [kother, base] = ipquad.base_other[i < k ? i : i - 1];
                kother = i;
                algoim::detail::eliminate_axis(ipquad.phi, kother, base.phi);
                base.build(false, true);
              }
          }
        }
      }

      // Makes phi exactly constant in the directions in which it is linear with (relative)
      // differences below tol (e.g. a plane that is parallel to a coordinate direction, up to
      // rounding errors). Algoim's root finding along such directions can produce NaNs otherwise.
      template <int N>
      void SnapConstantDirections (xarray<real,N> & alpha, real tol = 1.0e3 * std::numeric_limits<real>::epsilon())
      {
        real maxval = 0;
        for (int i = 0; i < alpha.size(); i++)
          maxval = std::max(maxval, std::abs(alpha[i]));
        for (int dim = 0; dim < N; dim++)
        {
          if (alpha.ext(dim) != 2)
            continue;
          real maxdiff = 0;
          for (auto i = alpha.loop(); ~i; ++i)
            if (i(dim) == 0)
            {
              uvector<int,N> j = i();
              j(dim) = 1;
              maxdiff = std::max(maxdiff, std::abs(alpha.m(j) - alpha.m(i())));
            }
          if (maxdiff > tol * maxval)
            continue;
          for (auto i = alpha.loop(); ~i; ++i)
            if (i(dim) == 0)
            {
              uvector<int,N> j = i();
              j(dim) = 1;
              const real mean = 0.5 * (alpha.m(j) + alpha.m(i()));
              alpha.m(i()) = mean;
              alpha.m(j) = mean;
            }
        }
      }

      // Adds the generic linear function eps * max|alpha| * (b + sum_d a_d y_d) to the Bernstein
      // polynomial alpha (Bernstein coefficients of a linear function are its values in the
      // Greville points). Several level sets in degenerate relative positions (e.g. planes that are
      // parallel to coordinate directions and intersection curves that coincide with each other or
      // with the boundary of the simplex) can lead to wrong quadrature rules in Algoim's elimination
      // steps; a tiny generic perturbation ("simulation of simplicity") avoids these degeneracies.
      template <int N>
      void PerturbGenerically (xarray<real,N> & alpha, int l, real eps)
      {
        // irrational coefficients (from the plastic number and the golden ratio), different for
        // every level set l
        const real a[3] = { 0.7548776662466927, 0.5698402909980532, 0.4301597090019468 };
        const real b = 0.3819660112501051;
        real maxval = 0;
        for (int i = 0; i < alpha.size(); i++)
          maxval = std::max(maxval, std::abs(alpha[i]));
        const real scale = eps * maxval;
        for (auto i = alpha.loop(); ~i; ++i)
        {
          real val = b * (l + 1);
          for (int d = 0; d < N; d++)
            val += a[(d + l) % 3] * (l + 1) * ((alpha.ext(d) > 1) ? real(i(d)) / (alpha.ext(d) - 1) : 0.5);
          alpha.m(i()) += scale * val;
        }
      }

      // sign condition of a level set: 0: phi < 0, 1: phi > 0, 2: phi = 0
      inline bool SignMatches (real val, int target)
      {
        return (target == 0) ? (val < 0) : (val > 0);
      }

      // Appends the rule for the part of B ∩ C where the level set functions phis[i] satisfy the
      // sign conditions targets[i] (at most one interface condition) to points/weights/gradients,
      // where B = lo + h [0,1]^N is a box in [0,1]^N and C is the cell (unit hypercube or unit
      // simplex). If Algoim detects vertical tangents in B (see above), B is subdivided into 2^N
      // boxes (at most 'levels' times).
      template <int N>
      void T_AppendBoxRule (const std::deque<xarray<real,N>> & phis, const std::vector<int> & targets,
                            bool simplex, const uvector<real,N> & lo, real h, int levels, int q,
                            std::vector<double> & points,
                            std::vector<double> & weights,
                            std::vector<double> & gradients)
      {
        using algoim::booluarray;
        const int n = phis.size();

        // the simplex {x_0 + ... + x_{N-1} < 1} is {psi < 0} with psi = y_0 + ... + y_{N-1} - L in
        // box coordinates y (x = lo + h y); the Bernstein coefficients of (multi-)linear psi are its
        // values in the corners of the unit cube. The mask restricts the quadrature algorithm to
        // the sub-cells (of Algoim's M x ... x M grid) that touch the closed simplex. Note that
        // Algoim assigns a root on the boundary of two sub-cells to the upper one (floor), so that
        // roots in a corner of a sub-cell are only taken into account if all sub-cells touching the
        // simplex are included. (lo and h are binary fractions, so all of this is exact.)
        const real L = (1 - sum(lo)) / h;
        xarray<real,N> psi(nullptr, uvector<int,N>(2));
        algoim_spark_alloc(real, psi);
        booluarray<N,ALGOIM_M> domain_mask(true);
        if (simplex)
        {
          for (MultiLoop<N> i(0, 2); ~i; ++i)
            psi.m(i()) = sum(i()) - L;
          for (MultiLoop<N> i(0, ALGOIM_M); ~i; ++i)
            domain_mask(i()) = sum(i()) <= L * ALGOIM_M;
        }
        auto inside = [&] (const uvector<real,N> & y) { return !simplex || (sum(y) < L); };

        // level set functions in box coordinates (exact, de Casteljau subdivision); level sets
        // without zero in B ∩ C either violate their sign condition (-> nothing to do) or are
        // satisfied in B ∩ C and need not be considered further
        const uvector<real,N> hi = lo + h;
        // (xarray is a non-copyable view, std::deque does not relocate its elements)
        std::vector<std::vector<real>> storage(n);
        std::deque<xarray<real,N>> phibs;
        std::vector<booluarray<N,ALGOIM_M>> masks;
        std::vector<int> cut_targets;
        int iface = -1; // index (in phibs) of the level set with interface condition
        const uvector<real,N> y0(real(0.5) / ALGOIM_M); // in B ∩ C (and in the first sub-cell)
        for (int i = 0; i < n; i++)
        {
          storage[i].resize(prod(phis[i].ext()));
          xarray<real,N> phib(storage[i].data(), phis[i].ext());
          algoim::bernstein::deCasteljau(phis[i], lo, hi, phib);
          if (n == 1)
            SnapConstantDirections<N>(phib);
          auto mask = algoim::detail::nonzeroMask(phib, domain_mask);
          if (algoim::detail::maskEmpty(mask))
          {
            if ((targets[i] == 2) || !SignMatches(algoim::bernstein::evalBernsteinPoly(phib, y0), targets[i]))
              return;
            continue;
          }
          if (targets[i] == 2)
            iface = phibs.size();
          phibs.emplace_back(storage[i].data(), phis[i].ext());
          masks.push_back(mask);
          cut_targets.push_back(targets[i]);
        }

        // quadrature hierarchy (dimension reduction w.r.t. suitable height directions); this is
        // what the constructors of ImplicitPolyQuadrature do, but with the masks computed above.
        // Without polynomials, a tensor-product Gauss rule is used.
        algoim::ImplicitPolyQuadrature<N> ipquad;
        for (size_t i = 0; i < phibs.size(); i++)
          ipquad.phi.push_back(phibs[i], masks[i]);
        if (simplex)
        {
          auto psi_mask = algoim::detail::nonzeroMask(psi, domain_mask);
          if (!algoim::detail::maskEmpty(psi_mask))
            ipquad.phi.push_back(psi, psi_mask);
        }
        if ((iface >= 0) && (phibs.size() > 1))
          BuildForInterface<N>(ipquad, phibs[iface], masks[iface]);
        else
          ipquad.build(true, false);

        if (!phibs.empty() && (levels > 0) && HasVerticalTangents<N>(ipquad))
        {
          for (MultiLoop<N> c(0, 2); ~c; ++c)
          {
            const uvector<real,N> lo_child = lo + (h/2) * uvector<real,N>(c());
            if (!simplex || (sum(lo_child) < 1))
              T_AppendBoxRule<N>(phis, targets, simplex, lo_child, h/2, levels-1, q, points, weights, gradients);
          }
          return;
        }

        // sign conditions of all (cut) level sets except for the interface condition
        auto signs_match = [&] (const uvector<real,N> & y)
        {
          for (size_t i = 0; i < phibs.size(); i++)
            if (((int)i != iface) && !SignMatches(algoim::bernstein::evalBernsteinPoly(phibs[i], y), cut_targets[i]))
              return false;
          return true;
        };

        if (iface >= 0)
        {
          if constexpr (N > 1)
          {
            // The surface rule contains points on the zero levels of all polynomials; the points on
            // {phi_iface = 0} are those for which this level set is the closest one (in terms of
            // |phi| / |grad phi|). Points on {psi = 0} (simplex) are discarded.
            const real wscale = std::pow(h, N-1);
            const real tol_face = 1e-12;
            auto dist = [] (const xarray<real,N> & p, const uvector<real,N> & y)
            {
              const real g = norm(algoim::bernstein::evalBernsteinPolyGradient(p, y));
              return std::abs(algoim::bernstein::evalBernsteinPoly(p, y)) / std::max(g, std::numeric_limits<real>::min());
            };
            ipquad.integrate_surf(strategy, q,
                                  [&] (const uvector<real,N> & y, real w, const uvector<real,N> & /*wn*/)
                                  {
                                    if ((w == 0) || (simplex && (sum(y) > L - tol_face)))
                                      return;
                                    if (phibs.size() > 1)
                                    {
                                      const real d = dist(phibs[iface], y);
                                      for (size_t i = 0; i < phibs.size(); i++)
                                        if (((int)i != iface) && (dist(phibs[i], y) < d))
                                          return;
                                      if (!signs_match(y))
                                        return;
                                    }
                                    // gradient w.r.t. y is h times the gradient w.r.t. x
                                    uvector<real,N> grad = algoim::bernstein::evalBernsteinPolyGradient(phibs[iface], y);
                                    for (int d = 0; d < N; d++)
                                    {
                                      points.push_back(lo(d) + h * y(d));
                                      gradients.push_back(grad(d));
                                    }
                                    weights.push_back(wscale * w);
                                  });
          }
          else
            throw std::logic_error("algoim_wrapper: 1D interfaces are treated separately");
        }
        else
        {
          const real wscale = std::pow(h, N);
          ipquad.integrate(strategy, q,
                           [&] (const uvector<real,N> & y, real w)
                           {
                             if ((w == 0) || !inside(y) || !signs_match(y))
                               return;
                             for (int d = 0; d < N; d++)
                               points.push_back(lo(d) + h * y(d));
                             weights.push_back(wscale * w);
                           });
        }
      }

      template <int N>
      MULTI_STATUS T_MultiCutQuadrature (int n, const int * degrees, const double * const * nodal_values,
                                         const int * targets_in, bool simplex, int q, int subdivlvl,
                                         std::vector<double> & points,
                                         std::vector<double> & weights,
                                         std::vector<double> & gradients)
      {
        using algoim::booluarray;
        points.clear();
        weights.clear();
        gradients.clear();

        booluarray<N,ALGOIM_M> domain_mask(true); // sub-cells of [0,1]^N that touch the cell
        if (simplex)
          for (MultiLoop<N> i(0, ALGOIM_M); ~i; ++i)
            domain_mask(i()) = sum(i()) <= ALGOIM_M;
        const uvector<real,N> center(simplex ? 1.0 / (N+1) : 0.5);

        std::vector<std::vector<real>> storage(n);
        std::deque<xarray<real,N>> phis; // views on storage (see T_AppendBoxRule)
        std::vector<int> targets;
        for (int l = 0; l < n; l++)
        {
          // Bernstein coefficients of the interpolating polynomial phi
          const int P = degrees[l] + 1;
          xarray<real,N> values(nullptr, uvector<int,N>(P));
          algoim_spark_alloc(real, values);
          storage[l].resize(prod(uvector<int,N>(P)));
          xarray<real,N> phi(storage[l].data(), uvector<int,N>(P));
          int k = 0;
          for (MultiLoop<N> i(0, P); ~i; ++i, ++k)
            values.m(i()) = nodal_values[l][k];
          algoim::bernstein::bernsteinInterpolate<N>(values, std::pow(100.0 * std::numeric_limits<real>::epsilon(), 1.0 / N), phi);
          // reduce the degree (per direction) if possible (up to a small tolerance), e.g. for a
          // quadratic level set function in a P3 space: the costs of Algoim's elimination steps
          // grow rapidly with the degree (resultants of degree O(degree^4) in 3D)
          ReduceDegree<N>(phi);
          if (n > 1)
            PerturbGenerically<N>(phi, l, multi_perturbation);
          else
            SnapConstantDirections<N>(phi);

          // level sets without zero in the cell: the sign condition holds everywhere or nowhere
          if (algoim::detail::maskEmpty(algoim::detail::nonzeroMask(phi, domain_mask)))
          {
            if ((targets_in[l] == 2) || !SignMatches(algoim::bernstein::evalBernsteinPoly(phi, center), targets_in[l]))
              return EMPTY;
            continue;
          }
          phis.emplace_back(storage[l].data(), phi.ext());
          targets.push_back(targets_in[l]);
        }
        if (phis.empty())
          return FULL;

        // one-dimensional interface: the roots of phi in (0,1) that satisfy the other conditions
        if constexpr (N == 1)
        {
          int iface = -1;
          for (size_t i = 0; i < phis.size(); i++)
            if (targets[i] == 2)
              iface = i;
          if (iface >= 0)
          {
            const int Pr = phis[iface].ext(0);
            std::vector<real> roots(std::max(Pr - 1, 1));
            const int nroots = algoim::bernstein::bernsteinUnitIntervalRealRoots(phis[iface].data(), Pr, roots.data());
            for (int r = 0; r < nroots; r++)
            {
              const uvector<real,1> x(roots[r]);
              bool ok = true;
              for (size_t i = 0; i < phis.size(); i++)
                if (((int)i != iface) && !SignMatches(algoim::bernstein::evalBernsteinPoly(phis[i], x), targets[i]))
                  ok = false;
              if (!ok)
                continue;
              points.push_back(roots[r]);
              weights.push_back(1.0);
              gradients.push_back(algoim::bernstein::evalBernsteinPolyGradient(phis[iface], x)(0));
            }
            return PARTIAL;
          }
        }

        // Uniform subdivision of [0,1]^N into nsub^N boxes (only those that intersect the cell are
        // used), each of them is adaptively subdivided once more if Algoim detects vertical
        // tangents (see HasVerticalTangents). Simplices are always subdivided (at least) once: this
        // keeps the boxes, in which Algoim works, comparable in size to the simplex, which
        // significantly improves the accuracy on coarse meshes (where a curved interface in a large
        // box tends to have (nearly) vertical tangents w.r.t. all directions). All subdivisions are
        // affine, hence the rules stay exact for polynomial integrands on planar cuts.
        const int nsub = (simplex ? 2 : 1) << subdivlvl;
        const int adaptive_levels = 1;
        for (MultiLoop<N> b(0, nsub); ~b; ++b)
          if (!simplex || (sum(b()) < nsub))
            T_AppendBoxRule<N>(phis, targets, simplex, (real(1) / nsub) * uvector<real,N>(b()), real(1) / nsub,
                               adaptive_levels, q, points, weights, gradients);
        return PARTIAL;
      }
    }

    void InterpolationPoints (int dim, int degree, std::vector<double> & points)
    {
      if (degree < 1)
        throw std::invalid_argument("algoim_wrapper::InterpolationPoints: degree has to be at least 1");
      switch (dim)
      {
      case 1: T_InterpolationPoints<1>(degree, points); break;
      case 2: T_InterpolationPoints<2>(degree, points); break;
      case 3: T_InterpolationPoints<3>(degree, points); break;
      default:
        throw std::invalid_argument("algoim_wrapper::InterpolationPoints: dim has to be 1, 2 or 3");
      }
    }

    MULTI_STATUS MultiCutQuadrature (int dim, int n, const int * degrees, const double * const * nodal_values,
                                     const int * domains, bool simplex, int q, int subdivlvl,
                                     std::vector<double> & points,
                                     std::vector<double> & weights,
                                     std::vector<double> & gradients)
    {
      int ninterfaces = 0;
      for (int l = 0; l < n; l++)
      {
        if (degrees[l] < 1)
          throw std::invalid_argument("algoim_wrapper::MultiCutQuadrature: degree has to be at least 1");
        if ((domains[l] < 0) || (domains[l] > 2))
          throw std::invalid_argument("algoim_wrapper::MultiCutQuadrature: domain has to be 0, 1 or 2");
        if (domains[l] == 2)
          ninterfaces++;
      }
      if (ninterfaces > 1)
        throw std::invalid_argument("algoim_wrapper::MultiCutQuadrature: at most one interface condition");
      if (q < 1)
        q = 1;
      if (subdivlvl < 0)
        subdivlvl = 0;
      switch (dim)
      {
      case 1: return T_MultiCutQuadrature<1>(n, degrees, nodal_values, domains, simplex, q, subdivlvl, points, weights, gradients);
      case 2: return T_MultiCutQuadrature<2>(n, degrees, nodal_values, domains, simplex, q, subdivlvl, points, weights, gradients);
      case 3: return T_MultiCutQuadrature<3>(n, degrees, nodal_values, domains, simplex, q, subdivlvl, points, weights, gradients);
      default:
        throw std::invalid_argument("algoim_wrapper::MultiCutQuadrature: dim has to be 1, 2 or 3");
      }
    }

    CUT_STATUS CutQuadrature (int dim, int degree, const double * nodal_values, bool simplex,
                              int domain, int q, int subdivlvl,
                              std::vector<double> & points,
                              std::vector<double> & weights,
                              std::vector<double> & gradients)
    {
      auto status = MultiCutQuadrature(dim, 1, &degree, &nodal_values, &domain, simplex, q, subdivlvl,
                                       points, weights, gradients);
      if (status == PARTIAL)
        return CUT;
      // EMPTY / FULL: the cell is on one side (for interfaces: an arbitrary side)
      const bool negative = (status == FULL) == (domain == 0);
      return negative ? UNCUT_NEGATIVE : UNCUT_POSITIVE;
    }
  }
}
