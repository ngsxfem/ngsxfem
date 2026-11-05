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

      // Appends the rule for the part 'domain' of B ∩ C to points/weights/gradients, where
      // B = lo + h [0,1]^N is a box in [0,1]^N and C is the cell (unit hypercube or unit simplex).
      // If Algoim detects vertical tangents in B (see above), B is subdivided into 2^N boxes
      // (at most 'levels' times).
      template <int N>
      void T_AppendBoxRule (const xarray<real,N> & phi, bool simplex,
                            const uvector<real,N> & lo, real h, int levels,
                            int domain, int q,
                            std::vector<double> & points,
                            std::vector<double> & weights,
                            std::vector<double> & gradients)
      {
        using algoim::booluarray;

        // level set function in box coordinates y, x = lo + h y (exact, de Casteljau subdivision)
        xarray<real,N> phib(nullptr, phi.ext());
        algoim_spark_alloc(real, phib);
        const uvector<real,N> hi = lo + h;
        algoim::bernstein::deCasteljau(phi, lo, hi, phib);

        // the simplex {x_0 + ... + x_{N-1} < 1} is {psi < 0} with psi = y_0 + ... + y_{N-1} - L in
        // box coordinates; the Bernstein coefficients of (multi-)linear psi are its values in the
        // corners of the unit cube. The mask restricts the quadrature algorithm to the sub-cells
        // (of Algoim's M x ... x M grid) that touch the closed simplex. Note that Algoim assigns a
        // root on the boundary of two sub-cells to the upper one (floor), so that roots in a
        // corner of a sub-cell are only taken into account if all sub-cells touching the simplex
        // are included. (lo and h are binary fractions, so all of this is exact.)
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

        auto phi_mask = algoim::detail::nonzeroMask(phib, domain_mask);
        const bool cut = !algoim::detail::maskEmpty(phi_mask);
        if (!cut)
        {
          // phi has no zero in B ∩ C: no interface, B ∩ C completely on one side
          if (domain == 2)
            return;
          uvector<real,N> y0(real(0.5) / ALGOIM_M); // in B ∩ C (and in the first sub-cell)
          if ((algoim::bernstein::evalBernsteinPoly(phib, y0) < 0) != (domain == 0))
            return;
        }

        // quadrature hierarchy (dimension reduction w.r.t. suitable height directions); this is
        // what the constructors of ImplicitPolyQuadrature do, but with the masks computed above.
        // Without polynomials, a tensor-product Gauss rule is used.
        algoim::ImplicitPolyQuadrature<N> ipquad;
        if (cut)
          ipquad.phi.push_back(phib, phi_mask);
        if (simplex)
        {
          auto psi_mask = algoim::detail::nonzeroMask(psi, domain_mask);
          if (!algoim::detail::maskEmpty(psi_mask))
            ipquad.phi.push_back(psi, psi_mask);
        }
        ipquad.build(true, false);

        if (cut && levels > 0 && HasVerticalTangents<N>(ipquad))
        {
          for (MultiLoop<N> c(0, 2); ~c; ++c)
          {
            const uvector<real,N> lo_child = lo + (h/2) * uvector<real,N>(c());
            if (!simplex || (sum(lo_child) < 1))
              T_AppendBoxRule<N>(phi, simplex, lo_child, h/2, levels-1, domain, q, points, weights, gradients);
          }
          return;
        }

        if (domain == 2)
        {
          if constexpr (N > 1)
          {
            // surface measure in x = surface measure in y times h^(N-1); for a simplex, the
            // surface rule also contains points on {psi = 0}, which are discarded
            const real wscale = std::pow(h, N-1);
            const real tol_face = 1e-12;
            ipquad.integrate_surf(strategy, q,
                                  [&] (const uvector<real,N> & y, real w, const uvector<real,N> & /*wn*/)
                                  {
                                    if ((w == 0) || (simplex && (sum(y) > L - tol_face)))
                                      return;
                                    // gradient w.r.t. y is h times the gradient w.r.t. x
                                    uvector<real,N> grad = algoim::bernstein::evalBernsteinPolyGradient(phib, y);
                                    for (int d = 0; d < N; d++)
                                    {
                                      points.push_back(lo(d) + h * y(d));
                                      gradients.push_back(grad(d));
                                    }
                                    weights.push_back(wscale * w);
                                  });
          }
          else
            throw std::invalid_argument("algoim_wrapper::CutQuadrature: no interface quadrature in 1D");
        }
        else
        {
          const real wscale = std::pow(h, N);
          const bool want_negative = (domain == 0);
          ipquad.integrate(strategy, q,
                           [&] (const uvector<real,N> & y, real w)
                           {
                             if ((w == 0) || !inside(y))
                               return;
                             const bool negative = algoim::bernstein::evalBernsteinPoly(phib, y) < 0;
                             if (negative != want_negative)
                               return;
                             for (int d = 0; d < N; d++)
                               points.push_back(lo(d) + h * y(d));
                             weights.push_back(wscale * w);
                           });
        }
      }

      template <int N>
      CUT_STATUS T_CutQuadrature (int degree, const double * nodal_values, bool simplex,
                                  int domain, int q, int subdivlvl,
                                  std::vector<double> & points,
                                  std::vector<double> & weights,
                                  std::vector<double> & gradients)
      {
        using algoim::booluarray;
        const int P = degree + 1;

        // Bernstein coefficients of the interpolating polynomial phi
        xarray<real,N> values(nullptr, uvector<int,N>(P));
        xarray<real,N> phi(nullptr, uvector<int,N>(P));
        algoim_spark_alloc(real, values, phi);
        int k = 0;
        for (MultiLoop<N> i(0, P); ~i; ++i, ++k)
          values.m(i()) = nodal_values[k];
        algoim::bernstein::bernsteinInterpolate<N>(values, std::pow(100.0 * std::numeric_limits<real>::epsilon(), 1.0 / N), phi);
        // reduce the degree (per direction) if possible (up to a small tolerance), e.g. for a
        // quadratic level set function in a P3 space: the costs of Algoim's elimination steps
        // grow rapidly with the degree (resultants of degree O(degree^4) in 3D)
        algoim::bernstein::autoReduction(phi);

        // check whether phi has a zero in the cell (sub-cells of [0,1]^N that touch the simplex)
        booluarray<N,ALGOIM_M> domain_mask(true);
        if (simplex)
          for (MultiLoop<N> i(0, ALGOIM_M); ~i; ++i)
            domain_mask(i()) = sum(i()) <= ALGOIM_M;
        if (algoim::detail::maskEmpty(algoim::detail::nonzeroMask(phi, domain_mask)))
        {
          uvector<real,N> center(simplex ? 1.0 / (N+1) : 0.5);
          return (algoim::bernstein::evalBernsteinPoly(phi, center) < 0) ? UNCUT_NEGATIVE : UNCUT_POSITIVE;
        }

        points.clear();
        weights.clear();
        gradients.clear();
        // Uniform subdivision of [0,1]^N into nsub^N boxes (only those that intersect the cell are
        // used), each of them is adaptively subdivided once more if Algoim detects vertical
        // tangents (see HasVerticalTangents). Simplices are always subdivided (at least) once: this keeps the boxes, in
        // which Algoim works, comparable in size to the simplex, which significantly improves the
        // accuracy on coarse meshes (where a curved interface in a large box tends to have
        // (nearly) vertical tangents w.r.t. all directions). All subdivisions are affine, hence
        // the rules stay exact for polynomial integrands on planar cuts.
        const int nsub = (simplex ? 2 : 1) << subdivlvl;
        const int adaptive_levels = 1;
        for (MultiLoop<N> b(0, nsub); ~b; ++b)
          if (!simplex || (sum(b()) < nsub))
            T_AppendBoxRule<N>(phi, simplex, (real(1) / nsub) * uvector<real,N>(b()), real(1) / nsub,
                               adaptive_levels, domain, q, points, weights, gradients);
        return CUT;
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

    CUT_STATUS CutQuadrature (int dim, int degree, const double * nodal_values, bool simplex,
                              int domain, int q, int subdivlvl,
                              std::vector<double> & points,
                              std::vector<double> & weights,
                              std::vector<double> & gradients)
    {
      if (degree < 1)
        throw std::invalid_argument("algoim_wrapper::CutQuadrature: degree has to be at least 1");
      if (q < 1)
        q = 1;
      if (subdivlvl < 0)
        subdivlvl = 0;
      switch (dim)
      {
      case 1: return T_CutQuadrature<1>(degree, nodal_values, simplex, domain, q, subdivlvl, points, weights, gradients);
      case 2: return T_CutQuadrature<2>(degree, nodal_values, simplex, domain, q, subdivlvl, points, weights, gradients);
      case 3: return T_CutQuadrature<3>(degree, nodal_values, simplex, domain, q, subdivlvl, points, weights, gradients);
      default:
        throw std::invalid_argument("algoim_wrapper::CutQuadrature: dim has to be 1, 2 or 3");
      }
    }
  }
}
