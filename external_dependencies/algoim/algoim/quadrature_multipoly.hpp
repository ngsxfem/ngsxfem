#ifndef ALGOIM_QUADRATURE_MULTIPOLY_HPP
#define ALGOIM_QUADRATURE_MULTIPOLY_HPP

// High-order accurate quadrature algorithms for multi-component domains implicitly
// defined by (one or more) multivariate Bernstein polynomials, based on the algorithms
// developed in the paper 
//    R. I. Saye, High-order quadrature on multi-component domains implicitly defined
//    by multivariate polynomials, Journal of Computational Physics, 448, 110720 (2022),
//    https://doi.org/10.1016/j.jcp.2021.110720
//
// See examples/examples_quad_multipoly.cpp, as well as the short tutorial on the GitHub
// page https://algoim.github.io/ for examples of usage.

// The algorithms make use of various "masking" operations. A mask divides the reference
// cube [0,1]^N into a regular grid of M x ... x M subcells; on each subcell, a mask has
// binary value 0 or 1 and indicates whether its accompanying polynomial is provably
// nonzero on that subcell, or if its roots can be ignored. Typically, M = 4 or 8 is
// a good choice, and this is specified by the following macro def.
#define ALGOIM_M 8

#include <algorithm>
#include "real.hpp"
#include "uvector.hpp"
#include "booluarray.hpp"
#include "multiloop.hpp"
#include "xarray.hpp"
#include "polyset.hpp"
#include "sparkstack.hpp"
#include "gaussquad.hpp"
#include "tanhsinh.hpp"
#include "bernstein.hpp"

namespace algoim
{
    namespace detail
    {
        template<int N>
        booluarray<N,ALGOIM_M> mask_driver(const xarray<real,N>& f, const booluarray<N,ALGOIM_M>& fmask, const xarray<real,N>* g, const booluarray<N,ALGOIM_M>* gmask)
        {
            booluarray<N,ALGOIM_M> mask(false);
            auto helper = [&](auto&& self, uvector<int,N> a, uvector<int,N> b)
            {
                bool overlap = false;
                for (MultiLoop<N> i(a,b); ~i; ++i)
                    if (fmask(i()) && (!gmask || (*gmask)(i())))
                        overlap = true;
                if (!overlap)
                    return;

                real eps = 0.015625 / ALGOIM_M;
                uvector<real,N> xa, xb;
                for (int dim = 0; dim < N; ++dim)
                {
                    xa(dim) = real(a(dim)) / ALGOIM_M - eps;
                    xb(dim) = real(b(dim)) / ALGOIM_M + eps;
                }

                if (g)
                {
                    xarray<real,N> fab(nullptr, f.ext()), gab(nullptr, g->ext());
                    algoim_spark_alloc(real, fab, gab);
                    bernstein::deCasteljau(f, xa, xb, fab);
                    bernstein::deCasteljau(*g, xa, xb, gab);
                    if (bernstein::orthantTest(fab, gab))
                        return;
                }
                else
                {
                    xarray<real,N> fab(nullptr, f.ext());
                    algoim_spark_alloc(real, fab);
                    bernstein::deCasteljau(f, xa, xb, fab);
                    if (bernstein::uniformSign(fab) != 0)
                        return;
                }

                if (b(0) - a(0) == 1)
                {
                    assert(all(b - a == 1));
                    assert(fmask(a) && (!gmask || (*gmask)(a)));
                    mask(a) = true;
                    return;
                }

                assert(all(b - a > 1) && all((b - a) % 2 == 0));
                uvector<int,N> delta = (b - a) / 2;
                for (MultiLoop<N> i(0,2); ~i; ++i)
                    self(self, a + i() * delta, a + (i() + 1) * delta);
            };
            helper(helper, 0, ALGOIM_M);
            return mask;
        }

        // Using orthant tests, compute a mask for the possible subrectangles for which f and g share
        // common zeros, i.e., intersecting zero level sets; if the mask is false somewhere, then it is
        // guaranteed that f and/or g were originally masked off at the same place, or that they
        // definitively do not share common zeros in that subrectangle; if the returned mask is true,
        // then shared zeros may exist (and with high likelihood)
        template<int N>
        booluarray<N,ALGOIM_M> intersectionMask(const xarray<real,N>& f, const booluarray<N,ALGOIM_M>& fmask, const xarray<real,N>& g, const booluarray<N,ALGOIM_M>& gmask)
        {
            return mask_driver(f, fmask, &g, &gmask);
        }

        // Using orthant tests, compute a mask for the possible subrectangles for which f has a zero
        // level set; if the mask is false somewhere, then it is guaranteed that f was originally
        // masked off at the same place, or that it definitively does not have any zeros in that
        // subrectangle; if the returned mask is true, then shared zeros may exist (and with high likelihood)
        template<int N>
        booluarray<N,ALGOIM_M> nonzeroMask(const xarray<real,N>& f, const booluarray<N,ALGOIM_M>& fmask)
        {
            return mask_driver<N>(f, fmask, nullptr, nullptr);
        }

        // Collapse a mask along dimension k by bitwise-OR-ing along columns
        template<int N>
        booluarray<N-1,ALGOIM_M> collapseMask(const booluarray<N,ALGOIM_M>& mask, int k)
        {
            booluarray<N-1,ALGOIM_M> r(false);
            for (MultiLoop<N> i(0,ALGOIM_M); ~i; ++i)
                if (mask(i()))
                    r(remove_component(i(), k)) = true;
            return r;
        }

        // Test if a mask is empty, i.e., all entries are false
        template<int N>
        bool maskEmpty(const booluarray<N,ALGOIM_M>& mask)
        {
            return mask.none();
        }

        // Test if a point x \in [0,1]^N is in a true subrectangle of a mask; if x is exactly
        // on the border between two subrectangles, the left subrectangle shall be used
        template<int N>
        bool pointWithinMask(const booluarray<N,ALGOIM_M>& mask, const uvector<real,N>& x)
        {
            using std::floor;
            uvector<int,N> cell;
            for (int dim = 0; dim < N; ++dim)
                cell(dim) = std::max(0, std::min(ALGOIM_M - 1, static_cast<int>(floor(x(dim) * ALGOIM_M))));
            return mask(cell);
        }

        // Test if a point {x + alpha e_k} is in a true subrectangle of a mask for some alpha \in [0,1]
        template<int N>
        bool lineIntersectsMask(const booluarray<N,ALGOIM_M>& mask, const uvector<real,N-1>& x, int k)
        {
            using std::floor;
            if constexpr (N > 1)
            {
                uvector<int,N> cell;
                for (int dim = 0; dim < N; ++dim)
                    if (dim < k)
                        cell(dim) = std::max(0, std::min(ALGOIM_M - 1, static_cast<int>(floor(x(dim) * ALGOIM_M))));
                    else if (dim > k)
                        cell(dim) = std::max(0, std::min(ALGOIM_M - 1, static_cast<int>(floor(x(dim - 1) * ALGOIM_M))));
                for (int i = 0; i < ALGOIM_M; ++i)
                {
                    cell(k) = i;
                    if (mask(cell))
                        return true;
                }
                return false;
            }
            else
                return !maskEmpty(mask);
        }

        template<int N>
        void restrictToFace(const xarray<real,N>& a, int k, int side, xarray<real,N-1>& out)
        {
            assert(0 <= k && k < N && (side == 0 || side == 1));
            assert(all(out.ext() == remove_component(a.ext(), k)));
            int P = a.ext(k);
            for (auto i = out.loop(); ~i; ++i)
            {
                uvector<int,N> j;
                for (int dim = 0; dim < N; ++dim)
                    j(dim) = (dim < k) ? i(dim) : ( (dim == k) ? side*(P-1) : i(dim - 1) );
                out.l(i) = a.m(j);
            }
        }

        template<int N>
        booluarray<N-1,ALGOIM_M> restrictToFace(const booluarray<N,ALGOIM_M>& a, int k, int side)
        {
            assert(0 <= k && k < N && (side == 0 || side == 1));
            booluarray<N-1,ALGOIM_M> r;
            for (MultiLoop<N-1> i(0,ALGOIM_M); ~i; ++i)
            {
                uvector<int,N> j;
                for (int dim = 0; dim < N; ++dim)
                    j(dim) = (dim < k)? i(dim) : ( (dim == k) ? side*(ALGOIM_M - 1) : i(dim - 1) );
                r(i()) = a(j);
            }
            return r;
        }

        // Compute determinant of the given matrix using QR + Givens rotations + column
        // pivoting, along with approximated rank
        //   in: square matrix A, which will be overwritten
        template<typename T>
        T det_qr(xarray<T,2>& A, int& rank, T tol = 10.0)
        {
            assert(A.ext(0) == A.ext(1) && A.ext(0) > 0);
            using std::max;
            using std::abs;
            T det = 1.0;
            int n = A.ext(0);
            T max_diag_r = 0.0;
            for (int j = 0; j < n; ++j)
            {
                T m = -1;
                int k = -1;
                for (int i = j; i < n; ++i)
                {
                    T mag = 0;
                    for (int a = 0; a < n; ++a)
                        mag += util::sqr(A(a,i));
                    if (k == -1 || mag >= m)
                    {
                        m = mag;
                        k = i;
                    }
                }
                assert(j <= k && k < n);
                if (k != j)
                {
                    for (int a = 0; a < n; ++a)
                        std::swap(A(a,j), A(a,k));
                    det *= -1.0;
                }
                for (int i = n - 1; i >= j + 1; --i)
                {
                    T c, s;
                    util::givens_get(A(i-1,j), A(i,j), c, s);
                    for (int k = j; k < n; ++k)
                        util::givens_rotate(A(i-1,k), A(i,k), c, s);
                }
                det *= A(j,j);
                max_diag_r = max(max_diag_r, abs(A(j,j)));
            }

            tol *= max_diag_r * n * std::numeric_limits<T>::epsilon();
            rank = 0;
            for (int i = 0; i < n; ++i)
                if (abs(A(i,i)) > tol)
                    ++rank;

            return det;
        }

        // Determine the largest possible degree of the resultant of two general polynomials
        template<int N>
        uvector<int,N-1> resultantExtent(const uvector<int,N>& p, const uvector<int,N>& q, int dim)
        {
            uvector<int,N-1> ext;
            for (int i = 0; i < N - 1; ++i)
            {
                int ii = (i < dim) ? i : i + 1;
                ext(i) = (p(dim) - 1) * (q(ii) - 1) + (q(dim) - 1) * (p(ii) - 1) + 1;
            }
            return ext;
        }

        // Determine the largest possible degree of the discriminant of a polynomial
        template<int N>
        uvector<int,N-1> discriminantExtent(const uvector<int,N>& p, int dim)
        {
            uvector<int,N-1> ext;
            for (int i = 0; i < N - 1; ++i)
            {
                int ii = (i < dim) ? i : i + 1;
                ext(i) = (2*p(dim) - 3) * (p(ii) - 1) + 1;
            }
            return ext;
        }

        // Compute the resultant of p and q and store the result in out
        // ========================================= NOTE =========================================
        // This is a heavily simplified method and does not handle rank deficiency caused by, e.g.,
        // common polynomial factors, nor does it handle ill-conditioning caused by extreme values,
        // among various other aspects. If your application requires this kind of special handling,
        // consider contacting the author of this code for suggestions.
        // ========================================= NOTE =========================================
        template<int N>
        bool resultant_core(const xarray<real,N>& p, const xarray<real,N>* q, int k, xarray<real,N-1>& out)
        {
            assert(0 <= k && k < N);

            int P = p.ext(k);
            int Q = q ? q->ext(k) : P - 1;
            int M = (P == Q) ? P - 1 : P + Q - 2;
            assert(P >= 2 && Q >= 1 && M >= 1);

            xarray<real,N-1> f(nullptr, out.ext());
            xarray<real,2> mat(nullptr, uvector<int,2>{M, M});
            real *pk, *qk;
            algoim_spark_alloc(real, f, mat);
            algoim_spark_alloc(real, &pk, P, &qk, Q);    

            for (auto i = f.loop(); ~i; ++i)
            {
                uvector<real,N-1> x;
                for (int dim = 0; dim < N - 1; ++dim)
                    x(dim) = bernstein::modifiedChebyshevNode(i(dim), f.ext(dim));
        
                bernstein::collapseAlongAxis(p, x, k, pk);
                if (q)
                    bernstein::collapseAlongAxis(*q, x, k, qk);
                else
                    bernstein::bernsteinDerivative(pk, P, qk);

                if (P == Q)
                    bernstein::bezoutMatrix(pk, qk, P, mat);
                else
                    bernstein::sylvesterMatrix(pk, P, qk, Q, mat);

                int rank;
                f.l(i) = det_qr(mat, rank);
            }

            // Interpolate the resultant on the tensor-product grid
            bernstein::normalise(f);
            bernstein::bernsteinInterpolate(f, std::pow(100.0 * std::numeric_limits<real>::epsilon(), 1.0 / (N - 1)), out);

            // Try for polynomial degree reduction
            bool b = bernstein::autoReduction(out, 1e4 * std::numeric_limits<real>::epsilon());

            // If able to reduce the degree, recompute the resultant on the lower-degree poly space
            // which is expected to have better conditioning
            if (b)
                resultant_core(p, q, k, out);
            return true;
        }

        // Compute the pseudo-resultant R(p,q) along dimension k
        template<int N>
        bool resultant(const xarray<real,N>& p, const xarray<real,N>& q, int k, xarray<real,N-1>& out)
        {
            return resultant_core(p, &q, k, out);
        }

        // Compute the (intentially unnormalised) pseudo-discriminant R(p,p') along dimension k
        template<int N>
        bool discriminant(const xarray<real,N>& p, int k, xarray<real,N-1>& out)
        {
            xarray<real,N> prime(nullptr, inc_component(p.ext(), k, -1));
            algoim_spark_alloc(real, prime);
            bernstein::bernsteinDerivative(p, k, prime);
            return resultant_core(p, &prime, k, out);
        }

        // Using the polynomials in phi, eliminate the axis k by restricting to faces, computing
        // discriminants and resultants, and storing the computed polynomials in psi
        template<int N>
        void eliminate_axis(PolySet<N,ALGOIM_M>& phi, int k, PolySet<N-1,ALGOIM_M>& psi)
        {
            static_assert(N >= 2, "N >= 2 required to eliminate axis");
            assert(0 <= k && k < N);
            assert(psi.count() == 0);

            // For every phi(i) ...
            for (int i = 0; i < phi.count(); ++i)
            {
                const auto& p = phi.poly(i);
                const auto& mask = phi.mask(i);

                // Examine bottom and top faces in the k'th dimension
                for (int side = 0; side <= 1; ++side)
                {
                    xarray<real,N-1> p_face(nullptr, remove_component(p.ext(), k));
                    algoim_spark_alloc(real, p_face);
                    restrictToFace(p, k, side, p_face);
                    auto p_face_mask = nonzeroMask(p_face, restrictToFace(mask, k, side));
                    if (!maskEmpty(p_face_mask))
                    {
                        bernstein::autoReduction(p_face);
                        bernstein::normalise(p_face);
                        psi.push_back(p_face, p_face_mask);
                    }
                }

                // Consider discriminant
                xarray<real,N> p_k(nullptr, p.ext());
                algoim_spark_alloc(real, p_k);
                bernstein::elevatedDerivative(p, k, p_k);
                auto disc_mask = intersectionMask(p, mask, p_k, mask);
                if (!maskEmpty(disc_mask))
                {
                    // note: computed disc might have lower degree than the following
                    uvector<int,N-1> R = discriminantExtent(p.ext(), k);
                    xarray<real,N-1> disc(nullptr, R);
                    algoim_spark_alloc(real, disc);
                    if (discriminant(p, k, disc))
                    {
                        bernstein::normalise(disc);
                        psi.push_back(disc, collapseMask(disc_mask, k));
                    }
                }
            }

            // Consider every pairwise combination of resultants ...
            for (int i = 0; i < phi.count(); ++i) for (int j = i + 1; j < phi.count(); ++j)
            {
                const auto& p = phi.poly(i);
                const auto& pmask = phi.mask(i);
                const auto& q = phi.poly(j);
                const auto& qmask = phi.mask(j);
                auto mask = intersectionMask(p, pmask, q, qmask);
                if (!maskEmpty(mask))
                {
                    // note: computed resultant might have lower degree than the following
                    uvector<int,N-1> R = resultantExtent(p.ext(), q.ext(), k);
                    xarray<real,N-1> res(nullptr, R);
                    algoim_spark_alloc(real, res);
                    if (resultant(p, q, k, res))
                    {
                        bernstein::normalise(res);
                        psi.push_back(res, collapseMask(mask, k));
                    }
                }
            };
        }

        // Compute the 'score' across all dimensions
        template<int N>
        uvector<real,N> score_estimate(PolySet<N,ALGOIM_M>& phi, uvector<bool,N>& has_disc)
        {
            static_assert(N > 1, "score_estimate of practical use only with N > 1");
            using std::abs;
            uvector<real,N> s = 0;
            has_disc = false;
            // For every phi(i) ...
            for (int i = 0; i < phi.count(); ++i)
            {
                const auto& p = phi.poly(i);
                const auto& mask = phi.mask(i);

                // Accumulate onto score by sampling at midpoint of every subcell of mask
                for (MultiLoop<N> j(0,ALGOIM_M); ~j; ++j) if (mask(j()))
                {
                    uvector<real,N> x = (j() + 0.5)/real(ALGOIM_M);
                    uvector<real,N> g = bernstein::evalBernsteinPolyGradient(p, x);
                    real sum = 0;
                    for (int dim = 0; dim < N; ++dim)
                    {
                        g(dim) = abs(g(dim));
                        sum += g(dim);
                    }
                    if (sum > 0)
                        s += g / sum;
                }

                // Consider discriminant
                xarray<real,N> p_k(nullptr, p.ext());
                algoim_spark_alloc(real, p_k);
                for (int k = 0; k < N; ++k)
                {
                    bernstein::elevatedDerivative(p, k, p_k);
                    auto disc_mask = intersectionMask(p, mask, p_k, mask);
                    has_disc(k) |= !maskEmpty(disc_mask);
                }
            }
            return s;
        }
    } // namespace detail

    // Main engine for generating high-order accurature quadrature schemes on multi-component domains
    //   implicitly defined by multivariate Bernstein polynomials in the unit hyperrectangle [0,1]^N.
    //   See examples_quad_multipoly.cpp for examples demonstrating the usage of these methods.
    //   See also additional examples provided on the GitHub documentation page,
    //   https://algoim.github.io/

    // After building the quadrature hierarchy (via dimension reduction), the specific kind
    // of quadrature scheme applied to the one-dimensional line integrals is chosen by the
    // user via a QuadStrategy parameter:
    //    AlwaysGL:  applies Gauss-Legendre quadrature on every interval, regardless of
    //               the level; this strategy is generally good for small-to-medium q
    //    AlwaysTS:  applies tanh-sinh quadrature on every interval, regardless of the level;
    //               this strategy is generally the worst of all three and is provided mainly
    //               for exploratory purposes
    //    AutoMixed: apply the automated strategy discussed in detail in the paper
    //               https://doi.org/10.1016/j.jcp.2021.110720; briefly, the method: (i) applies
    //               tanh-sinh on base integrals for which vertical tangents are detected in
    //               the height-function-based representation of their implicitly-defined
    //               geometry; (ii) applies Gauss-Legendre quadrature on the inner-most
    //               integral, and on all base integrals whose geometry is the graph of a 
    //               multi-valued height function devoid of vertical tangents/branching.
    enum QuadStrategy { AlwaysGL, AlwaysTS, AutoMixed };

    template<int N>
    struct ImplicitPolyQuadrature
    {
        enum IntegralType { Inner, OuterSingle, OuterAggregate };

        PolySet<N,ALGOIM_M> phi;                                                // Given N-dimensional polynomials
        int k;                                                                  // Elimination axis/height direction; k = N if there are no interfaces
        ImplicitPolyQuadrature<N-1> base;                                       // Base polynomials corresponding to removal of axis k
        bool auto_apply_TS;                                                     // If quad method is auto chosen, indicates whether TS is applied
        IntegralType type;                                                      // Whether an inner integral, or outer of two kinds
        std::array<std::tuple<int,ImplicitPolyQuadrature<N-1>>,N-1> base_other; // Stores other base cases, besides k, when in aggregate mode

        // Default ctor sets to an uninitialised state
        ImplicitPolyQuadrature() : k(-1) {}

        // Build quadrature hierarchy for a domain implicitly defined by a single polynomial
        ImplicitPolyQuadrature(const xarray<real,N>& p)
        {
            auto mask = detail::nonzeroMask(p, booluarray<N,ALGOIM_M>(true));
            if (!detail::maskEmpty(mask))
                phi.push_back(p, mask);
            build(true, false);
        }

        // Build quadrature hierarchy for a domain implicitly defined by two polynomials
        ImplicitPolyQuadrature(const xarray<real,N>& p, const xarray<real,N>& q)
        {
            {
                auto mask = detail::nonzeroMask(p, booluarray<N,ALGOIM_M>(true));
                if (!detail::maskEmpty(mask))
                    phi.push_back(p, mask);
            }
            {
                auto mask = detail::nonzeroMask(q, booluarray<N,ALGOIM_M>(true));
                if (!detail::maskEmpty(mask))
                    phi.push_back(q, mask);
            }
            build(true, false);
        }

        // Build quadrature hierarchy for a given domain implicitly defined by two polynomials with user-defined masks
        ImplicitPolyQuadrature(const xarray<real,N>& p, const booluarray<N,ALGOIM_M>& pmask, const xarray<real,N>& q, const booluarray<N,ALGOIM_M>& qmask)
        {
            {
                auto mask = detail::nonzeroMask(p, pmask);
                if (!maskEmpty(mask))
                    phi.push_back(p, mask);
            }
            {
                auto mask = detail::nonzeroMask(q, qmask);
                if (!maskEmpty(mask))
                    phi.push_back(q, mask);
            }
            build(true, false);
        }

        // Assuming phi has been instantiated, determine elimination axis and build base
        void build(bool outer, bool auto_apply_TS)
        {
            type = outer ? OuterSingle : Inner;
            this->auto_apply_TS = auto_apply_TS;

            // If phi is empty, apply a tensor-product Gaussian quadrature
            if (phi.count() == 0)
            {
                k = N;
                this->auto_apply_TS = false;
                return;
            }

            if constexpr (N == 1)
            {
                // If in one dimension, there is only one choice of height direction and
                // the recursive process halts
                k = 0;
                return;
            }
            else
            {
                // Compute score; penalise any directions which likely contain vertical tangents
                uvector<bool,N> has_disc;
                uvector<real,N> score = detail::score_estimate(phi, has_disc);
                //if (max(abs(score)) == 0)
                //    score(0) = 1.0;
                assert(max(abs(score)) > 0);
                score /= 2 * max(abs(score));
                for (int i = 0; i < N; ++i)
                    if (!has_disc(i))
                        score(i) += 1.0;

                // Choose height direction and form base polynomials; if tanh-sinh is being used at this
                // level, suggest the same all the way down; moreover, suggest tanh-sinh if a non-empty
                // discriminant mask has been found
                k = argmax(score);
                detail::eliminate_axis(phi, k, base.phi);
                base.build(false, this->auto_apply_TS || has_disc(k));

                // If this is the outer integral, and surface quadrature schemes are required, apply
                // the dimension-aggregated scheme when necessary
                if (outer && has_disc(k))
                {
                    type = OuterAggregate;
                    for (int i = 0; i < N; ++i) if (i != k)
                    {
                        auto& [kother, base] = base_other[i < k ? i : i - 1];
                        kother = i;
                        detail::eliminate_axis(phi, kother, base.phi);
                        // In aggregate mode, triggered by non-empty discriminant mask,
                        // base integrals always have T-S suggested
                        base.build(false, this->auto_apply_TS || true);
                    }
                }
            }
        }

        // Integrate a functional via quadrature of the base integral, adding the dimension k
        template<typename F>
        void integrate(QuadStrategy strategy, int q, const F& f)
        {
            assert(0 <= k && k <= N);

            // If there are no interfaces, apply rectangular tensor-product Gauss-Legendre quadrature
            if (k == N)
            {
                assert(!auto_apply_TS);
                for (MultiLoop<N> i(0, q); ~i; ++i)
                {
                    uvector<real,N> x;
                    real w = 1.0;
                    for (int dim = 0; dim < N; ++dim)
                    {
                        x(dim) = GaussQuad::x(q, i(dim));
                        w *= GaussQuad::w(q, i(dim));
                    }
                    f(x, w);
                }
                return;
            }

            // Determine maximum possible number of roots; used to allocate a small buffer
            int max_count = 2;
            for (size_t i = 0; i < phi.count(); ++i)
                max_count += phi.poly(i).ext(k) - 1;

            // Base integral invokes the following integrand
            auto integrand = [&](const uvector<real,N-1>& xbase, real w)
            {
                // Allocate node buffer of sufficient size and initialise with {0, 1}
                real *nodes;
                algoim_spark_alloc(real, &nodes, max_count);
                nodes[0] = 0.0;
                nodes[1] = 1.0;
                int count = 2;

                // For every phi(i) ...
                for (size_t i = 0; i < phi.count(); ++i)
                {
                    const auto& p = phi.poly(i);
                    const auto& mask = phi.mask(i);
                    int P = p.ext(k);

                    // Ignore phi if its mask is void everywhere above the base point
                    if (!detail::lineIntersectsMask(mask, xbase, k))
                        continue;

                    // Restrict polynomial to axis-aligned line and compute its roots
                    real *pline, *roots;
                    algoim_spark_alloc(real, &pline, P, &roots, P - 1);
                    bernstein::collapseAlongAxis(p, xbase, k, pline);
                    int rcount = bernstein::bernsteinUnitIntervalRealRoots(pline, P, roots);

                    // Add all real roots in [0,1] which are also within masked region of phi
                    for (int j = 0; j < rcount; ++j)
                    {
                        auto x = add_component(xbase, k, roots[j]);
                        if (detail::pointWithinMask(mask, x))
                            nodes[count++] = roots[j];
                    }
                };

                // Sort the nodes
                std::sort(nodes, nodes + count);
                assert(nodes[0] == real(0) && nodes[count-1] == real(1));

                // Force nearly-degenerate sub-intervals to be exactly degenerate
                real tol = 10.0 * std::numeric_limits<real>::epsilon();
                using std::abs;
                for (int i = 1; i < count - 1; ++i)
                    if (abs(nodes[i]) < tol)
                        nodes[i] = 0.0;
                    else if (abs(nodes[i] - 1) < tol)
                        nodes[i] = 1.0;
                    else if (abs(nodes[i] - nodes[i+1]) < tol)
                        nodes[i+1] = nodes[i];
                assert(nodes[0] == real(0) && nodes[count-1] == real(1));

                // Apply quadrature to non-degenerate sub-intervals
                for (int i = 0; i < count - 1; ++i)
                {
                    real x0 = nodes[i];
                    real x1 = nodes[i + 1];
                    if (x0 == x1)
                        continue;

                    // Choose between Gauss-Legendre and tanh-sinh according to the user-defined strategy
                    bool GL = true;
                    if (strategy == AlwaysTS)
                        GL = false;
                    if (strategy == AutoMixed)
                        GL = !this->auto_apply_TS;
                
                    if (GL)
                        for (int j = 0; j < q; ++j)
                            f(add_component(xbase, k, x0 + (x1 - x0) * GaussQuad::x(q, j)), w * (x1 - x0) * GaussQuad::w(q, j));
                    else
                        for (int j = 0; j < q; ++j)
                            f(add_component(xbase, k, TanhSinhQuadrature::x(q, j, x0, x1)), w * TanhSinhQuadrature::w(q, j, x0, x1));
                }
            };

            // When N = 1, the base case of recursion is invoked on a zero-dimensional point with unit weight
            if constexpr (N > 1)
                base.integrate(strategy, q, integrand);
            else
                integrand(uvector<real,0>(), real(1));
        }

        // Surface-integrate a functional via quadrature of the base integral, adding the dimension k
        template<typename F>
        void integrate_surf(QuadStrategy strategy, int q, const F& f)
        {
            static_assert(N > 1, "surface integral only implemented in N > 1 dimensions");
            assert(type == OuterSingle || type == OuterAggregate);

            // If there is no interface, there is no surface integral
            if (k == N)
                return;

            // Base integral invokes the following integrand which operates in the height direction k_active
            int k_active = -1;
            auto integrand = [&](const uvector<real,N-1>& xbase, real w)
            {
                assert(0 <= k_active && k_active < N);
                // For every phi(i) ...
                for (size_t i = 0; i < phi.count(); ++i)
                {
                    const auto& p = phi.poly(i);
                    const auto& mask = phi.mask(i);
                    int P = p.ext(k_active);

                    // Ignore phi if its mask is void everywhere above the base point
                    if (!detail::lineIntersectsMask(mask, xbase, k_active))
                        continue;

                    // Compute roots of { x \mapsto phi(xbase + x e_k) }
                    real *pline, *roots;
                    algoim_spark_alloc(real, &pline, P, &roots, P - 1);
                    bernstein::collapseAlongAxis(p, xbase, k_active, pline);
                    int rcount = bernstein::bernsteinUnitIntervalRealRoots(pline, P, roots);

                    // Consider all real roots in (0,1) which are also within masked region of phi; evaluate
                    // integrand at interfacial points, multiplying weights by the effective surface Jacobian
                    for (int j = 0; j < rcount; ++j)
                    {
                        auto x = add_component(xbase, k_active, roots[j]);
                        if (detail::pointWithinMask(mask, x))
                        {
                            using std::abs;
                            uvector<real,N> g = bernstein::evalBernsteinPolyGradient(p, x);
                            if (type == OuterAggregate)
                            {
                                // When in aggregate mode, the scalar surf integral multiplies f by |n_k|^2, whose net effect
                                // is multiply weight by |n_k|; the flux surf integral multiplies f by sign(n_k) = sign(g(k))
                                real alpha = max(abs(g));
                                if (alpha > 0)
                                {
                                    g /= alpha;
                                    alpha = abs(g(k_active)) / norm(g);
                                }
                                // Simplistic method to compute sign(n_k). NOTE: This method relies on a reasonable consistency
                                // between the gradient calculation of original polynomial, and that of the roots computed from
                                // pline; when near high-multiplicity roots, this simple method can break down; other, more
                                // sophisticated methods should be used in such cases, but these are not implemented here
                                f(x, w * alpha, set_component<real,N>(real(0.0), k_active, w * util::sign(g(k_active))));
                            }
                            else
                            {
                                // When in non-aggregate mode, the scalar surf integral multiples f by 1, whose net effect
                                // is multiply weight by 1/|n_k|; the flux surf integral multiplies f by n
                                uvector<real,N> n = g;
                                if (norm(n) > 0)
                                    n *= real(1.0) / norm(n);
                                real alpha = w * norm(g) / abs(g(k_active));
                                f(x, alpha, alpha * n);
                            }
                        }
                    }
                };
            };

            // Apply primary base integral
            k_active = k;
            base.integrate(strategy, q, integrand);

            // If in aggregate mode, apply to other dimensions as well
            if (type == OuterAggregate)
            {
                for (int i = 0; i < N - 1; ++i)
                {
                    auto& [k, base] = base_other[i];
                    k_active = k;
                    base.integrate(strategy, q, integrand);
                }
            }
        }
    };

    template<> struct ImplicitPolyQuadrature<0> {};
} // namespace algoim

#endif
