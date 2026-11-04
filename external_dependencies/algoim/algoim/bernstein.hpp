#ifndef ALGOIM_BERNSTEIN_HPP
#define ALGOIM_BERNSTEIN_HPP

// algoim::bernstein implements several routines for working with multivariate Bernstein
// polynomials. Many of these methods, especially the orthant, root finding, Sylvester
// and Bezout methods, are based on those described in the paper
//    R. I. Saye, High-order quadrature on multi-component domains implicitly defined
//    by multivariate polynomials, Journal of Computational Physics, 448, 110720 (2022),
//    https://doi.org/10.1016/j.jcp.2021.110720

#include <cassert>
#include "real.hpp"
#include "uvector.hpp"
#include "xarray.hpp"
#include "sparkstack.hpp"
#include "binomial.hpp"
#include "utility.hpp"

// Some methods rely on a LAPACK implementation to solve
// generalised eigenvalue problems and SVD factorisation
#if __has_include(<lapacke.h>)
#include <lapacke.h>
#elif __has_include(<mkl_lapacke.h>)
#include <mkl_lapacke.h>
#else
#error "Algoim requires a LAPACKE implementation to compute eigenvalues and SVD factorisations, but a suitable lapacke.h include file was not found; did you forget to specify its include path?"
#endif

namespace algoim::bernstein
{
    // Evaluate at x the P Bernstein basis functions of degree P-1
    //   out: array of length P
    template<typename T>
    void evalBernsteinBasis(const T& x, int P, T* out)
    {
        assert(P >= 1);
        const real* binom = Binomial::row(P - 1);
        T p = 1.0;
        for (int i = 0; i < P; ++i)
        {
            out[i] = p * binom[i];
            p *= x;
        }
        p = 1.0;
        for (int i = P - 1; i >= 0; --i)
        {
            out[i] *= p;
            p *= 1.0 - x;
        }
    }

    // Evaluate an N-dimensional Bernstein polynomial at x
    template<int N>
    real evalBernsteinPoly(const xarray<real,N>& beta, const uvector<real,N>& x)
    {
        uvector<real*,N> basis;
        algoim_spark_alloc_vec(real, basis, beta.ext());
        for (int i = 0; i < N; ++i)
            evalBernsteinBasis(x(i), beta.ext(i), basis(i));
        real r = 0.0;
        for (auto i = beta.loop(); ~i; ++i)
        {
            real s = beta.l(i);
            for (int dim = 0; dim < N; ++dim)
                s *= basis(dim)[i(dim)];
            r += s;
        }
        return r;
    }

    // Fast evaluation of a 1-D Bernstein polynomial and its derivative; it is assumed that
    // binom == Binomial::row(P - 1), left to the caller to evaluate and cache, for speed
    void bernsteinValueAndDerivative(const real* alpha, int P, const real* binom, real x, real& value, real& deriv)
    {
        assert(P > 1);
        real *a, *b;
        algoim_spark_alloc(real, &a, P, &b, P);
        a[0] = 1; for (int i = 1; i < P; ++i) a[i] = a[i-1] * x;
        b[0] = 1; for (int i = 1; i < P; ++i) b[i] = b[i-1] * (1 - x);
        value = alpha[0] * b[P-1] + alpha[P-1] * a[P-1];
        for (int i = 1; i < P - 1; ++i)
            value += alpha[i] * binom[i] * a[i] * b[P-1-i];
        deriv = (alpha[P-1] * a[P-2] - alpha[0] * b[P-2]) * (P - 1);
        for (int i = 1; i < P-1; ++i)
            deriv += alpha[i] * binom[i] * (a[i-1] * b[P-1-i] * i - a[i] * b[P-2-i] * (P-1-i));
    }

    // Evaluate the gradient of an N-dimensional Bernstein polynomial at x
    template<int N>
    uvector<real,N> evalBernsteinPolyGradient(const xarray<real,N>& beta, const uvector<real,N>& x)
    {
        uvector<real*,N> basis, prime;
        algoim_spark_alloc_vec(real, basis, beta.ext());
        algoim_spark_alloc_vec(real, prime, beta.ext());
        for (int i = 0; i < N; ++i)
        {
            int P = beta.ext(i);
            assert(P >= 1);
            evalBernsteinBasis(x(i), P, basis(i));
            if (P > 1)
            {
                real *buff;
                algoim_spark_alloc(real, &buff, P - 1);
                evalBernsteinBasis(x(i), P - 1, buff);
                prime(i)[0]     = (P - 1) * (          - buff[0]);
                prime(i)[P-1]   = (P - 1) * (buff[P-2]          );
                for (int j = 1; j < P - 1; ++j)
                    prime(i)[j] = (P - 1) * (buff[j-1] - buff[j]);
            }
            else
                prime(i)[0] = 0.0;
        }
        uvector<real,N> g = real(0.0);
        for (auto i = beta.loop(); ~i; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                real s = beta.l(i);
                for (int dim = 0; dim < N; ++dim)
                    if (dim == j)
                        s *= prime(dim)[i(dim)];
                    else
                        s *= basis(dim)[i(dim)];
                g(j) += s;
            }
        }
        return g;
    }

    // Assuming p is represented in scaled Bernstein coefficients, reverse that scaling
    template<int N>
    void reverseScaledCoeff(xarray<real,N>& p)
    {
        uvector<const real*,N> binom;
        for (int i = 0; i < N; ++i)
            binom(i) = Binomial::row(p.ext(i) - 1);
        for (auto i = p.loop(); i; ++i)
        {
            real alpha = 1;
            for (int dim = 0; dim < N; ++dim)
                alpha *= binom(dim)[i(dim)];
            p.l(i) /= alpha;
        }
    }

    // Squared L2 norm of a Bernstein polynomial; the result may be negative, but only if
    // the polynomial is essentially machine zero
    template<int N>
    real squaredL2norm(const xarray<real,N>& p)
    {
        uvector<const real*,N> b1, b2;
        for (int dim = 0; dim < N; ++dim)
        {
            b1(dim) = Binomial::row(p.ext(dim) - 1);
            b2(dim) = Binomial::row(2*p.ext(dim) - 2);
        }
        real delta = 0;
        for (auto i = p.loop(); ~i; ++i) for (auto j = p.loop(); ~j; ++j)
        {
            real g = 1;
            for (int dim = 0; dim < N; ++dim)
                g *= (b1(dim)[i(dim)] / b2(dim)[i(dim) + j(dim)]) * b1(dim)[j(dim)];
            delta += p.l(i) * p.l(j) * g;
        }
        for (int dim = 0; dim < N; ++dim)
            delta /= 2*p.ext(dim) - 1;
        return delta;
    }

    // Collapse a multivariate Bernstein polynomial along a given axis-aligned line, i.e., the
    // set of points x where x(dim) is free and x(i) = x0(i) for all i != dim
    //   out: array of length beta.ext(dim)
    template<int N>
    void collapseAlongAxis(const xarray<real,N>& beta, const uvector<real,N-1>& x0, int dim, real* out)
    {
        if constexpr (N == 1)
        {
            assert(dim == 0);
            for (int i = 0; i < beta.ext(0); ++i)
                out[i] = beta[i];
        }
        else
        {
            assert(0 <= dim && dim < N);
            uvector<real*,N-1> basis;
            algoim_spark_alloc_vec(real, basis, remove_component(beta.ext(), dim));
            for (int i = 0; i < N - 1; ++i)
            {
                int P = beta.ext(i < dim ? i : i + 1);
                evalBernsteinBasis(x0(i), P, basis(i));
            }
            int P = beta.ext(dim);
            for (int i = 0; i < P; ++i)
                out[i] = 0.0;
            for (auto i = beta.loop(); ~i; ++i)
            {
                real s = beta.l(i);
                for (int j = 0; j < N; ++j)
                    if (j < dim)
                        s *= basis(j)[i(j)];
                    else if (j > dim)
                        s *= basis(j-1)[i(j)];
                out[i(dim)] += s;
            }
        }
    }

    // Collapse a multivariate Bernstein polynomial along a given axis-orthogonal hyperplane,
    // i.e., the set of points x where x(dim) is a fixed given value
    template<int N>
    void collapseAlongHyperplane(const xarray<real,N>& beta, int dim, real x, xarray<real,N-1>& out)
    {
        static_assert(N > 1, "N > 1 required");
        assert(all(out.ext() == remove_component(beta.ext(), dim)));
        assert(0 <= dim && dim < N);
        int P = beta.ext(dim);
        real *basis;
        algoim_spark_alloc(real, &basis, P);
        evalBernsteinBasis(x, P, basis);
        out = 0.0;
        for (auto i = beta.loop(); i; ++i)
            out.m(remove_component(i(), dim)) += beta.l(i) * basis[i(dim)];
    }

    // Normalise polynomial by its largest (in absolute value) coefficient
    template<int N>
    void normalise(xarray<real,N>& alpha)
    {
        real x = alpha.maxNorm();
        if (x > 0)
            alpha *= real(1) / x;
    }

    // Applying a simple examination of coefficient signs, returns +1 if the
    // polynomial is uniformly positive, -1 if the polynomial is uniformly
    // negative, or 0 if no guarantees can be made
    template<int N>
    int uniformSign(const xarray<real,N>& beta) 
    {
        int s = util::sign(beta[0]);
        for (int i = 1; i < beta.size(); ++i)
            if (util::sign(beta[i]) != s)
                return 0;
        return s;
    }

    // Compute coefficients of derivative in lower degree basis
    //   alpha: array of length P
    //   out: array of length P - 1
    template<typename T>
    void bernsteinDerivative(const T* alpha, int P, T* out)
    {
        assert(P >= 2);
        for (int i = 0; i < P - 1; ++i)
        {
            out[i] = alpha[i + 1];
            out[i] -= alpha[i];
            out[i] *= P - 1;
        }
    }

    // Compute the derivative of a Bernstein polynomial
    template<int N>
    void bernsteinDerivative(const xarray<real,N>& a, int dim, xarray<real,N>& out)
    {
        assert(all(out.ext() == inc_component(a.ext(), dim, -1)));
        int P = a.ext(dim);
        assert(P >= 2);
        for (auto i = out.loop(); ~i; ++i)
            out.l(i) = a.m(i.shifted(dim, 1)) - a.m(i());
        out *= P - 1;
    }

    // Compute the derivative of a Bernstein polynomial, in the original basis;
    // equivalent to computing normal derivative, and then elevating once
    template<int N>
    void elevatedDerivative(const xarray<real,N>& a, int dim, xarray<real,N>& out)
    {
        assert(all(out.ext() == a.ext()) && 0 <= dim && dim < N);
        int P = a.ext(dim);
        for (auto i = a.loop(); ~i; ++i)
        {
            if (i(dim) == 0)
                out.l(i) = (a.m(i.shifted(dim, 1)) - a.l(i)) * (P - 1);
            else if (i(dim) == P - 1)
                out.l(i) = (a.l(i) - a.m(i.shifted(dim, -1))) * (P - 1);
            else
                out.l(i) = a.m(i.shifted(dim, -1)) * (-i(dim)) + a.l(i) * (2*i(dim) - P + 1) + a.m(i.shifted(dim, 1)) * (P - 1 - i(dim));
        }
    }

    // Apply de Casteljau algorithm to compute the Bernstein coefficients of alpha relative to the interval [0,tau]
    template<int N>
    void deCasteljauLeft(xarray<real,N>& alpha, real tau)
    {
        int P = alpha.ext(0);
        for (int i = 1; i < P; ++i)
            for (int j = P - 1; j >= i; --j)
            {
                alpha.a(j) *= tau;
                alpha.a(j) += alpha.a(j - 1) * (1.0 - tau);
            }
    }

    // Apply de Casteljau algorithm to compute the Bernstein coefficients of alpha relative to the interval [tau,1]
    template<int N>
    void deCasteljauRight(xarray<real,N>& alpha, real tau)
    {
        int P = alpha.ext(0);
        for (int i = 1; i < P; ++i)
            for (int j = 0; j < P - i; ++j)
            {
                alpha.a(j) *= (1.0 - tau);
                alpha.a(j) += alpha.a(j + 1) * tau;
            }
    }

    // Apply de Casteljau algorithm to compute the Bernstein coefficients of alpha relative
    // to the hyperrectangle [a,b]. It is assumed that the given arrays a & b each have
    // length at least N. If, for a particular dimension, a[dim] > b[dim], both the interval
    // and coefficients are reversed.
    template<int N, bool B = false>
    void deCasteljau(xarray<real,N>& alpha, const real* a, const real* b)
    {
        using std::swap;
        using std::abs;
        if constexpr (N == 1 || B)
        {
            int P = alpha.ext(0);
            if (*b < *a)
            {
                deCasteljau<N,B>(alpha, b, a);
                for (int i = 0; i < P / 2; ++i)
                    swap(alpha.a(i), alpha.a(P - 1 - i));
                return;
            }
            if (abs(*b) >= abs(*a - 1))
            {
                deCasteljauLeft(alpha, *b);
                deCasteljauRight(alpha, *a / *b);
            }
            else
            {
                deCasteljauRight(alpha, *a);
                deCasteljauLeft(alpha, (*b - *a)/(real(1) - *a));
            }
        }
        else
        {
            deCasteljau<2,true>(alpha.flatten().ref(), a, b);
            for (int i = 0; i < alpha.ext(0); ++i)
                deCasteljau(alpha.slice(i).ref(), a + 1, b + 1);
        }
    }

    // Apply de Casteljau algorithm to compute the Bernstein coefficients of alpha relative
    // to the hyperrectangle [a,b]. If, for a particular dimension, a[dim] > b[dim], both the
    // interval and coefficients are reversed.
    template<int N>
    void deCasteljau(const xarray<real,N>& alpha, const uvector<real,N>& a, const uvector<real,N>& b, xarray<real,N>& out)
    {
        assert(all(out.ext() == alpha.ext()));
        out = alpha;
        deCasteljau(out, a.data(), b.data());
    }

    // Elevate the degree of a Bernstein polynomial
    template<int N, bool B = false>
    void bernsteinElevate(const xarray<real,N>& alpha, xarray<real,N>& beta)
    {
        assert(all(beta.ext() >= alpha.ext()));
        if constexpr (N == 1 || B)
        {
            int P = alpha.ext(0), Q = beta.ext(0);
            if (P == Q)
            {
                for (int k = 0; k < P; ++k)
                    beta.a(k) = alpha.a(k);
            }
            else
            {
                int n = P - 1;
                int r = Q - 1 - n;
                if (r == 1)
                {
                    beta.a(0) = alpha.a(0);
                    beta.a(n+1) = alpha.a(n);
                    for (int k = 1; k <= n; ++k)
                    {
                        beta.a(k) = alpha.a(k - 1) * (real(k) / real(n + 1));
                        beta.a(k) += alpha.a(k) * (real(1) - real(k) / real(n + 1));
                    }
                    return;
                }
                const real* bn = Binomial::row(n);
                const real* br = Binomial::row(r);
                const real* bnr = Binomial::row(n + r);
                for (int k = 0; k <= n + r; ++k)
                {
                    beta.a(k) = 0.0;
                    for (int j = std::max(0, k - r); j <= std::min(n, k); ++j)
                        beta.a(k) += alpha.a(j) * ((br[k-j] * bn[j]) / bnr[k]);
                }
            }
        }
        else
        {
            xarray<real,N> gamma(nullptr, set_component(alpha.ext(), 0, beta.ext(0)));
            algoim_spark_alloc(real, gamma);
            bernsteinElevate<2,true>(alpha.flatten(), gamma.flatten().ref());
            for (int i = 0; i < beta.ext(0); ++i)
                bernsteinElevate(gamma.slice(i), beta.slice(i).ref());
        }
    }

    namespace detail
    {
        // Compute least squares solution of Ax=b, where A is a (P+1) x P lower bidiagonal matrix, with
        // diagonal given by 'alpha', lower diagonal by 'beta', each of length P, and the rhs is given
        // by 'b', of length P+1; the algorithm applies QR with Givens, which shall create an upper
        // bidiagonal R; no pivoting is performed, and the complexity is O(P) for QR factorisation and
        // O(P*O) for back-solve on rhs of size PxO
        //   in: alpha, length P, shall be overwritten
        //   in: beta, length P, shall be overwritten
        //   in: b, overwritten with first P rows yielding least squares solution
        void lsqr_bidiagonal(real *alpha, real *beta, int P, xarray<real,2>& b)
        {
            assert(b.ext(0) == P + 1 && b.ext(1) > 0);
            real *gamma;
            algoim_spark_alloc_def(real, 0, &gamma, P);
            for (int i = 0; i < P; ++i)
            {
                real c, s;
                util::givens_get(alpha[i], beta[i], c, s);
                util::givens_rotate(alpha[i], beta[i], c, s);
                if (i < P - 1) util::givens_rotate(gamma[i+1], alpha[i+1], c, s);
                for (int k = 0; k < b.ext(1); ++k)
                    util::givens_rotate(b(i,k), b(i+1,k), c, s);
            }
            b.a(P-1) /= alpha[P-1];
            for (int i = P - 2; i >= 0; --i)
            {
                b.a(i) -= b.a(i + 1) * gamma[i + 1];
                b.a(i) /= alpha[i];
            }
        }
    } // namespace detail

    // Reduce by one the effective degree of a Bernstein polynomial; this routine
    // mainly makes sense when the actual polynomial degree is less than the one
    // used in its (starting) Bernstein polynomial representation
    template<int N, bool B = false>
    void bernsteinReduction(xarray<real,N>& alpha, int dim)
    {
        assert(all(alpha.ext() >= 1) && 0 <= dim && dim < N && alpha.ext(dim) >= 2);
        if (dim == 0)
        {
            int P = alpha.ext(0) - 1;
            real *a, *b;
            algoim_spark_alloc(real, &a, P, &b, P);
            a[0] = 1;
            b[P-1] = 1;
            for (int k = 1; k < P; ++k)
            {
                a[k] = real(1) - real(k) / real(P);
                b[k-1] = real(k) / real(P);
            }
            xarray<real,2> view(alpha.data(), uvector<int,2>{P + 1, prod(alpha.ext(), 0)});
            detail::lsqr_bidiagonal(a, b, P, view);
        }
        else if constexpr (N > 1)
        {
            for (int i = 0; i < alpha.ext(0); ++i)
                bernsteinReduction<N-1,true>(alpha.slice(i).ref(), dim - 1);
        }

        if (!B)
        {
            xarray<real,N> beta(nullptr, alpha.ext());
            algoim_spark_alloc(real, beta);
            beta = alpha;
            alpha.alterExtent(inc_component(alpha.ext(), dim, -1));
            for (auto i = alpha.loop(); ~i; ++i)
                alpha.l(i) = beta.m(i());
        }
    }

    // Automatically reduce the degree of alpha; returns true iff degree reduction occurred
    template<int N>
    bool autoReduction(xarray<real,N>& alpha, real tol = 1.0e3 * std::numeric_limits<real>::epsilon(), int dim = 0)
    {
        using std::sqrt;
        using std::abs;
        if (dim < 0 || dim >= N || tol <= 0)
            return false;
        bool stay = false;
        if (alpha.ext(dim) >= 2)
        {
            xarray<real,N> beta(nullptr, alpha.ext()), gamma(nullptr, alpha.ext());
            algoim_spark_alloc(real, beta, gamma);
            beta = alpha;
            bernsteinReduction(beta, dim);
            bernsteinElevate(beta, gamma);
            gamma -= alpha;
            real delta = sqrt(abs(squaredL2norm(gamma)));
            real norm = sqrt(abs(squaredL2norm(alpha)));
            if (delta < tol * norm)
            {
                alpha.alterExtent(beta.ext());
                alpha = beta;
                stay = true;
            }
        }
        if (stay)
        {
            autoReduction<N>(alpha, tol, dim);
            return true;
        }
        else
            return autoReduction<N>(alpha, tol, dim + 1);
    }

    // Determine if there is a scalar alpha such that sign x(i) + alpha y(i) > 0 for every component i;
    // if sign = 0, then returns true if it holds for sign = 1 and/or sign = -1
    template<int N>
    bool orthantTestBase(const xarray<real,N>& x, const xarray<real,N>& y, int sign = 0)
    {
        assert(sign == 0 || sign == -1 || sign == 1);
        assert(all(x.ext() == y.ext()));
        using std::min;
        using std::max;
        using std::isinf;
        using std::abs;
        if (sign == 0)
            return orthantTestBase(x, y, -1) || orthantTestBase(x, y, 1);
        real alphaMax =  std::numeric_limits<real>::infinity();
        real alphaMin = -std::numeric_limits<real>::infinity();
        for (int i = 0; i < x.size(); ++i)
        {
            if (y[i] == 0.0 && x[i] * sign <= 0.0)
                return false;
            if (y[i] > 0.0)
                alphaMin = max(alphaMin, -x[i] / y[i] * sign);
            else if (y[i] < 0.0)
                alphaMax = min(alphaMax, -x[i] / y[i] * sign);
        }
        if (isinf(alphaMin) || isinf(alphaMax))
            return true;
        if (alphaMax - alphaMin > 1.0e5 * std::numeric_limits<real>::epsilon() * max(abs(alphaMin), abs(alphaMax)))
            return true;
        return false;
    }

    // Determine if there are scalars alpha and beta such that {alpha f + beta g > 0} holds for every
    // Bernstein coefficient of f and g: if one of the polynomials has a smaller degree than the other,
    // it is degree elevated so that the two polynomials have the same degree
    template<int N>
    bool orthantTest(const xarray<real,N>& f, const xarray<real,N>& g)
    {
        if (all(f.ext() == g.ext()))
            return orthantTestBase(f, g);
        else
        {
            uvector<int,N> ext = max(f.ext(), g.ext());
            xarray<real,N> fe(nullptr, ext), ge(nullptr, ext);
            algoim_spark_alloc(real, fe, ge);
            bernsteinElevate(f, fe);
            bernsteinElevate(g, ge);
            return orthantTestBase(fe, ge);
        }
    }

    // Modified Chebyshev nodes (which include endpoints) for interpolating degree P-1 polynomials
    inline real modifiedChebyshevNode(int i, int P)
    {
        assert(0 <= i && i < P);
        using std::cos;
        return 0.5 - 0.5 * cos(util::pi * i / (P - 1));
    }

    // Methods to compute, and cache, the SVD for Bernstein interpolation based on modified Chebysev nodes
    struct BernsteinVandermondeSVD
    {
        struct SVD
        {
            real *U;
            real *Vt;
            real *sigma;
        };

        static SVD get(int P)
        {
            assert(P >= 1);
            static thread_local std::unordered_map<int,std::vector<real>> cache;
            if (cache.count(P) == 1)
            {
                real *base = cache.at(P).data();
                return SVD{base, base + P*P, base + 2*P*P};
            }

            real *A, *superb, *basis;
            algoim_spark_alloc(real, &A, P*P, &superb, P, &basis, P);
            for (int i = 0; i < P; ++i)
            {
                evalBernsteinBasis(modifiedChebyshevNode(i, P), P, basis);
                for (int j = 0; j < P; ++j)
                    A[i*P + j] = basis[j];
            }

            cache[P].resize(P*P + P*P + P);
            real *base = cache[P].data();
            SVD result{base, base + P*P, base + 2*P*P};

            static_assert(std::is_same_v<real, double>, "Algoim's default LAPACK code assumes real == double; a custom SVD solver is required when real != double");
            int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', P, P, A, P, result.sigma, result.U, P, result.Vt, P, superb);
            assert(info == 0 && "LAPACKE_dgesvd call failed (algoim::bernstein::BernsteinVandermondeSVD::get)");
            return result;
        }
    };

    // Interpolate tensor-product data f, assumed to be nodal values at the same nodes returned by modifiedChebyshevNode()
    template<int N, bool B = false>
    void bernsteinInterpolate(const xarray<real,N>& f, real tol, xarray<real,N>& out)
    {
        assert(all(out.ext() == f.ext()));
        if constexpr (N == 1 || B)
        {
            int P = f.ext(0);
            int O = prod(f.ext(), 0);
            assert(P >= 1 && O >= 1);

            real *tmp;
            algoim_spark_alloc(real, &tmp, P * O);

            auto svd = BernsteinVandermondeSVD::get(P);

            for (int i = 0; i < P * O; ++i)
                tmp[i] = 0.0;
            for (int i = 0; i < P; ++i)
                for (int j = 0; j < P; ++j)
                    for (int k = 0; k < O; ++k)
                        tmp[i*O + k] += svd.U[j*P + i] * f[j*O + k];

            real minsigma = tol * svd.sigma[0];
            for (int i = 0; i < P; ++i)
            {
                real alpha = (svd.sigma[i] >= minsigma) ? (real(1) / svd.sigma[i]) : 0.0;
                for (int k = 0; k < O; ++k)
                    tmp[i*O + k] *= alpha;
            }

            out = 0;
            for (int i = 0; i < P; ++i)
                for (int j = 0; j < P; ++j)
                    for (int k = 0; k < O; ++k)
                        out[i*O + k] += svd.Vt[j*P + i] * tmp[j*O + k];
        }
        else
        {
            xarray<real,N> gamma(nullptr, f.ext());
            algoim_spark_alloc(real, gamma);        
            bernsteinInterpolate<2,true>(f.flatten(), tol, gamma.flatten().ref());
            for (int i = 0; i < f.ext(0); ++i)
                bernsteinInterpolate(gamma.slice(i), tol, out.slice(i).ref());
        }
    }

    // Interpolate a functional through its nodal evaluation at the modifiedChebyshevNode() points
    template<int N, typename F>
    void bernsteinInterpolate(F&& f, xarray<real,N>& out)
    {
        xarray<real,N> ff(nullptr, out.ext());
        algoim_spark_alloc(real, ff);
        for (auto i = ff.loop(); ~i; ++i)
        {
            uvector<real,N> x;
            for (int dim = 0; dim < N; ++dim)
                x(dim) = modifiedChebyshevNode(i(dim), out.ext(dim));
            ff.l(i) = f(x);
        }
        bernsteinInterpolate(ff, std::pow(100.0 * std::numeric_limits<real>::epsilon(), 1.0 / N), out);
    }

    namespace detail
    {
        // Compute the generalised eigenvalues for matrix pair A, B
        //   in: N by N square matrices; A, B will be overwritten
        //   out: array of length N x 2
        void generalisedEigenvalues(xarray<real,2>& A, xarray<real,2>& B, xarray<real,2>& out)
        {
            int N = A.ext(0);
            assert(all(A.ext() == N) && all(B.ext() == N) && out.ext(0) == N && out.ext(1) == 2);
            real *alphar, *alphai, *beta, *lscale, *rscale;
            algoim_spark_alloc(real, &alphar, N, &alphai, N, &beta, N, &lscale, N, &rscale, N);
            real abnrm, bbnrm;
            int ilo, ihi;
            static_assert(std::is_same_v<real, double>, "Algoim's default LAPACK code assumes real == double; a custom generalised eigenvalue solver is required when real != double");
            int info = LAPACKE_dggevx(LAPACK_ROW_MAJOR, 'B', 'N', 'N', 'N', N, A.data(), N, B.data(), N, alphar, alphai, beta, nullptr, N, nullptr, N, &ilo, &ihi, lscale, rscale, &abnrm, &bbnrm, nullptr, nullptr);
            assert(info == 0 && "LAPACKE_dggevx call failed (algoim::bernstein::detail::generalisedEigenvalues)");
            for (int i = 0; i < N; ++i)
            {
                if (beta[i] != 0.0)
                    out(i,0) = alphar[i] / beta[i],
                    out(i,1) = alphai[i] / beta[i];
                else
                    out(i,0) = out(i,1) = std::numeric_limits<real>::infinity();
            }
        }
    } // namespace detail

    // Compute all complex  roots of a Bernstein polynomial
    //   alpha: array of length P
    //   out: array of length (P-1) x 2
    void rootsBernsteinPoly(const real* alpha, int P, xarray<real,2>& out)
    {
        assert(P >= 2 && out.ext(0) == P - 1 && out.ext(1) == 2);
        using std::max;
        using std::abs;

        real *beta;
        algoim_spark_alloc(real, &beta, P);
        real tol = 0.0;
        for (int i = 0; i < P; ++i)
            tol = max(tol, abs(alpha[i]));
        tol *= util::sqr(std::numeric_limits<real>::epsilon());
        for (int i = 0; i < P; ++i)
            beta[i] = (abs(alpha[i]) > tol) ? alpha[i] : 0;

        int N = P - 1;
        xarray<real,2> A(nullptr, uvector<int,2>{N, N});
        xarray<real,2> B(nullptr, uvector<int,2>{N, N});
        algoim_spark_alloc(real, A, B);
        A = 0;
        B = 0;
        for (int i = 0; i < N - 1; ++i)
            A(i, i + 1) = B(i, i + 1) = 1.0;
        for (int i = 0; i < N; ++i)
            A(N - 1, i) = B(N - 1, i) = -beta[i];
        B(N - 1, N - 1) += beta[N] / N;
        for (int i = 0; i < N - 1; ++i)
            B(i, i) = real(N - i) / real(i + 1);

        detail::generalisedEigenvalues(A, B, out);
    }

    namespace detail
    {
        // Newton's method safeguarded by a standard bisection method; in Bernstein application it
        // is only be applied to a Bernstein polynomial guaranteed to have just one real root
        template<typename F>
        bool newtonBisectionSearch(const F& f, real x0, real x1, real tol, int maxsteps, real& root)
        {
            using std::abs;
            real f0, f1, dummy;
            f(x0, f0, dummy);
            f(x1, f1, dummy);
            if ((f0 > 0.0 && f1 > 0.0) || (f0 < 0.0 && f1 < 0.0))
                return false;
            if (f0 == real(0.0))
            {
                root = x0;
                return true;
            }
            if (f1 == real(0.0))
            {
                root = x1;
                return true;
            }

            // x0 and x1 define the bracket; x0 always corresponds to negative value of f; x1 positive value of f
            if (f1 < 0.0)
                std::swap(x0, x1);

            // Initial guess is midpoint
            real x = (x0 + x1)*0.5;
            real fx, fpx;
            f(x, fx, fpx);
            real dx = x1 - x0;
            for (int step = 0; step < maxsteps; ++step)
            {
                if ((fpx*(x - x0) - fx)*(fpx*(x - x1) - fx) < 0.0 && abs(fx) < abs(dx*fpx)*0.5)
                {
                    // Step in Newton's method falls within bracket and is less than half the previous step size
                    dx = -fx / fpx;
                    real xold = x;
                    x += dx;
                    if (xold == x)
                    {
                        root = x;
                        return true;
                    }
                }
                else
                {
                    // Revert to bisection
                    dx = (x1 - x0)*0.5;
                    x = x0 + dx;
                    if (x == x0)
                    {
                        root = x;
                        return true;
                    }
                }
                if (abs(dx) < tol)
                {
                    root = x;
                    return true;
                }
                f(x, fx, fpx);
                if (fx == real(0.0)) // Got very lucky
                {
                    root = x;
                    return true;
                }
                if (fx < 0.0)
                    x0 = x;
                else
                    x1 = x;
            }
            return false;
        }
    }

    // Compute, if possible, a simple real root in [0,1] of a Bernstein polynomial using 
    // Descartes' rule of signs:
    //   - if it can be guaranteed that there is exactly 0 roots, 0 is returned
    //   - if it can be guaranteed that there is exactly 1 root, and that root has
    //     been calculated to full precision using Newton's method, 1 is returned
    //   - if some coefficients are close to zero (thereby preventing a reliable use of
    //     Descartes' rule), -1 is returned
    //   - if no other guarantees can be made, -1 is returned
    int bernsteinSimpleRoot(const real* alpha, int P, real tol, real& root)
    {
        assert(P >= 2);
        using std::abs;
        for (int i = 0; i < P; ++i)
            if (abs(alpha[i]) < tol)
                return -1;
        int count = 0;
        for (int i = 1; i < P; ++i)
            if (alpha[i-1] < 0 && alpha[i] >= 0 || alpha[i-1] >= 0 && alpha[i] < 0)
                ++count;
        if (count == 0)
            return 0;
        if (count > 1)
            return -1;
        real newton_tol = 10.0 * std::numeric_limits<real>::epsilon();
        const real* binom = Binomial::row(P - 1);
        bool b = detail::newtonBisectionSearch([=](real x, real& value, real& prime)
        {
            bernsteinValueAndDerivative(alpha, P, binom, x, value, prime);
        }, 0, 1, newton_tol, 12, root);
        return b ? 1 : -1;
    }

    // Compute real roots of a Bernstein polynomial using a bisection + Newton's method approach.
    // Returns the number of real roots computed (and recorded in out, a buffer of size at least
    // P - 1), or -1 if failed
    int rootsBernsteinPolyFast(const xarray<real,1>& alpha, real a, real b, int depth, real tol, real* out)
    {
        // Try simple root method
        real root;
        int res = bernsteinSimpleRoot(alpha.data(), alpha.ext(0), tol, root);
        // If it worked with a guarantee of no roots, return
        if (res == 0)
            return 0;
        // If it worked with a guarantee of just one root computed accurately,
        // transform that root to the [a,b] interval, record it, and return
        if (res == 1)
        {
            *out = a + (b - a)*root;
            return 1;
        }
        // Otherwise, the simple root method failed. Apply bisection, provided not already too deep
        if (depth >= 4)
            return -1;
        xarray<real,1> beta(nullptr, alpha.ext());
        algoim_spark_alloc(real, beta);
        // Apply to left half
        beta = alpha;
        deCasteljauLeft(beta, 0.5);
        int r1 = rootsBernsteinPolyFast(beta, a, a + (b - a)*0.5, depth + 1, tol, out);
        if (r1 < 0)
            return -1;
        // Apply to right half, shifting buffer by r1 
        beta = alpha;
        deCasteljauRight(beta, 0.5);
        int r2 = rootsBernsteinPolyFast(beta, a + (b - a)*0.5, b, depth + 1, tol, out + r1);
        if (r2 < 0)
            return -1;
        return r1 + r2;
    }

    // Apply generalised eigenvalue method to compute the real roots of alpha in the interval [0,1],
    // returning the number of roots recorded in 'out', a buffer of size at least P - 1
    int bernsteinUnitIntervalRealRoots_eigenvalue(const real* alpha, int P, real* out)
    {
        using std::abs;
        xarray<real,2> roots(nullptr, uvector<int,2>{P - 1, 2});
        algoim_spark_alloc(real, roots);
        rootsBernsteinPoly(alpha, P, roots);
        real tol = 1.0e4 * std::numeric_limits<real>::epsilon(); // nearly-real-root tolerance
        int count = 0;
        for (int j = 0; j < P - 1; ++j)
        {
            if (0 <= roots(j,0) && roots(j,0) <= 1 && abs(roots(j,1)) < tol)
            {
                *(out + count) = roots(j,0);
                ++count;
            }
        }
        return count;
    }

    // Apply a Newton's method-based approach to compute the real roots of alpha in the interval [0,1];
    // if succeeded, returns the number of roots recorded in 'out' (a buffer of size at least P -1);
    // if failed, returns -1
    int bernsteinUnitIntervalRealRoots_fast(const real* alpha, int P, real* out)
    {
        using std::max;
        using std::abs;
        // Compute a tolerance by which to declare a nearly-zero coefficient as being
        // too close to zero (for Descartes' rule of signs and to avoid problems where
        // a root lies close to a subinterval endpoint which can confuse bisection)
        real tol = 0;
        for (int i = 0; i < P; ++i)
            tol = max(tol, abs(alpha[i]));
        tol *= 1.0e4 * std::numeric_limits<real>::epsilon(); // nearly zero coeff tolerance, can be loose
        return rootsBernsteinPolyFast(xarray<real,1>(const_cast<real*>(alpha), P), 0, 1, 0, tol, out);
    }

    // Driver method to compute the real roots of a Bernstein polynomial in the interval [0,1];
    // the method first tries a fast approach, which succeeds in the vast majority of cases and is
    // anywhere between 10x and 100x faster than the backup approach; if the fast approach fails,
    // the backup method is applied. Returns the number of computed roots, recorded in the buffer
    // 'out' of size at least P - 1
    int bernsteinUnitIntervalRealRoots(const real* alpha, int P, real* out)
    {
        using std::sqrt;
        if (P == 1)
            return 0;

        // Direct method for linear polynomials
        if (P == 2)
        {
            if (alpha[0] == alpha[1]) return 0;
            real x = alpha[0] / (alpha[0] - alpha[1]);
            if (x < 0 || x > 1) return 0;
            *out = x;
            return 1;
        }

        // Direct method for quadratic polynomials, using numerically-stable quadratic formula
        if (P == 3)
        {
            real a = alpha[0] - alpha[1] * 2 + alpha[2];
            real b = (alpha[1] - alpha[0]) * 2;
            real c = alpha[0];
            real delta = b*b - a*c*4;
            if (delta < 0) return 0;
            real q = -0.5 * (b + (b >= 0 ? sqrt(delta) : -sqrt(delta)));
            real r1 = q / a;
            real r2 = c / q;
            int count = 0;
            if (0 <= r1 && r1 <= 1) { *out = r1; ++count; }
            if (0 <= r2 && r2 <= 1) { *(out + count) = r2; ++count; }
            return count;
        }

        // Apply fast method, if possible, and resort to eigenvalue method if it fails
        int count = bernsteinUnitIntervalRealRoots_fast(alpha, P, out);
        if (count >= 0)
            return count;
        return bernsteinUnitIntervalRealRoots_eigenvalue(alpha, P, out);
    }

    // Build Sylvester matrix for Bernstein polynomials of degrees P-1 and Q-1
    //   out: square matrix of dimensions P + Q - 2
    void sylvesterMatrix(const real* a, int P, const real* b, int Q, xarray<real,2>& out)
    {
        assert(P >= 1 && Q >= 1 && P + Q >= 3 && out.ext(0) == P + Q - 2 && out.ext(1) == P + Q - 2);
        const real* bP = Binomial::row(P - 1);
        const real* bQ = Binomial::row(Q - 1);
        const real* bPQ = Binomial::row(P + Q - 3);
        out = 0;
        for (int i = 0; i < Q - 1; ++i)
            for (int j = 0; j < P; ++j)
                out(i,j+i) = a[j] * (bP[j] / bPQ[j + i]);
        for (int i = 0; i < P - 1; ++i)
            for (int j = 0; j < Q; ++j)
                out(i+Q-1,j+i) = b[j] * (bQ[j] / bPQ[j + i]);
    }

    // Build Bezout matrix for Bernstein polynomials of equal degree P-1
    //   out: square matrix of dimensions P - 1
    void bezoutMatrix(const real* a, const real* b, int P, xarray<real,2>& out)
    {
        assert(P >= 2 && out.ext(0) == P - 1 && out.ext(1) == P - 1);
        const int n = P - 1;
        out = 0;
        for (int i = 1; i <= n; ++i)
            out(i-1,0) = (a[i] * b[0] - a[0] * b[i]) * real(n) / real(i);
        for (int j = 1; j <= n - 1; ++j)
            out(n-1,j) = (a[n] * b[j] - a[j] * b[n]) * real(n) / real(n - j);
        for (int i = n - 1; i >= 1; --i)
            for (int j = 1; j <= i - 1; ++j)
                out(i-1,j) = (a[i] * b[j] - a[j] * b[i]) * real(n*n) / real(i*(n - j)) + out(i,j-1) * real(j*(n-i)) / real(i*(n-j));
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                out(i,j) = out(j,i);
    }
} // namespace algoim::bernstein

#endif
