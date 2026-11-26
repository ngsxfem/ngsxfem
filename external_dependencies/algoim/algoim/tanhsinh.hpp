#ifndef ALGOIM_TANHSINH_HPP
#define ALGOIM_TANHSINH_HPP

// algoim::TanhSinhQuadrature implements basic tanh-sinh quadrature methods. Some
// of the design choices follow the discussion in the paper
//    R. I. Saye, High-order quadrature on multi-component domains implicitly defined
//    by multivariate polynomials, Journal of Computational Physics, 448, 110720 (2022),
//    https://doi.org/10.1016/j.jcp.2021.110720

#include <cassert>
#include <array>
#include <vector>
#include "real.hpp"
#include "uvector.hpp"

namespace algoim
{
    struct TanhSinhQuadrature
    {
        static constexpr int n_max = 100;

        // Computes and caches all quadrature points and nodes for schemes
        // containing up to n_max many quadrature points
        static const std::array<real,n_max*(n_max+1)>& data()
        {
            auto generator = []()
            {
                std::array<real,n_max*(n_max+1)> d;
                for (int n = 1; ; ++n)
                {
                    // Generate a scheme using at most n points
                    std::vector<real> scheme(2*n);
                    int m = generate(n, &scheme[0]);
                    // m is the number of effective nodes; store the scheme
                    // which has the largest n for an effective m
                    if (m > n_max + 10)
                        break;
                    if (m <= n_max)
                        for (int i = 0; i < 2*m; ++i)
                            d[m*(m-1) + i] = scheme[i];
                }
                return d;
            };
            static auto data_ = generator();
            return data_;
        };

        // Quadrature node, relative to the interval [-1,1], for an n-point scheme
        static real x(int n, int i)
        {
            assert(1 <= n && n <= n_max && 0 <= i && i < n);
            return data()[n*(n-1) + 2*i + 0];
        }

        // Quadrature weight, relative to the interval [-1,1], for an n-point scheme
        static real w(int n, int i)
        {
            assert(1 <= n && n <= n_max && 0 <= i && i < n);
            return data()[n*(n-1) + 2*i + 1];
        }

        // Quadrature node, relative to the interval [a,b], for an n-point scheme
        static real x(int n, int i, real a, real b)
        {
            return (a + b + (b - a) * x(n, i)) * 0.5;
        }

        // Quadrature weight, relative to the interval [a,b], for an n-point scheme
        static real w(int n, int i, real a, real b)
        {
            return w(n, i) * (b - a) * 0.5;
        }

        // Generate quadrature points and weights for an n-point scheme; the result is stored
        // in pairs (x,w) in the output buffer, required to be of length at least 2n. Depending
        // on a tolerance, one or more endpoint nodes may be snapped to -1 or 1, and if more
        // than one, their quadrature weights accumulated; the result returned is the
        // reduced number of effective nodes
        static int generate(int n, real* out)
        {
            assert(n >= 1);

            // midpoint rule for one quadrature point
            if (n == 1)
            {
                out[0] = 0.0;
                out[1] = 2;
                return 1;
            }

            // Newton's method to compute the Lambert W function W(z), accurate to full precision
            auto LambertW = [](real z)
            {
                using std::log;
                using std::exp;
                real w = z < 1.0 ? z - 0.45 * z * z : 0.75 * log(z);
                // 6 iterations needed for double; 8 generally enough for quad-double;
                // definitely possible to improve with a more sophisticated initial guess
                for (int i = 0; i < 10; ++i)
                    w -= (w * exp(w) - z) / (exp(w) + w * exp(w));
                return w;
            };

            real h_a = 0.6 * 0.5 * util::pi;
            real h = 2.0 * LambertW(2.0 * h_a * (n - 1)) / n;

            int count = 0;

            // Generate central node when number of points is odd
            if (n % 2)
            {
                out[count++] = 0.0;
                out[count++] = util::pi * 0.5;
            }

            bool snappedEndpoint = false;

            // Generate remaining points symmetrically in pairs
            using std::exp;
            using std::abs;
            for (int i = 0; i < n/2; ++i)
            {
                real t = (n % 2) ? (i + 1) * h : (i + 0.5) * h;
                real exp_t = exp(t);
                real exp_mt = 1.0 / exp_t;
                real omega = 0.25 * util::pi * (exp_t - exp_mt); // omega = 0.5 pi sinh(t)
                real exp_omega = exp(omega);
                real cosh_omega_2 = exp_omega + 1.0 / exp_omega; // 2 cosh(omega)
                real cosh_t_2 = exp_t + exp_mt; // 2 cosh(t)

                // The quadrature weight function w(t) satisfies
                //   w(t) = pi (2 cosh(t)) / (2 cosh(omega))^2
                // The complementary quadrature node, 1 - x, satisfies
                //   1 - x(t) = 2 / (1 + exp(2 omega))
                real w = util::pi * cosh_t_2 / (cosh_omega_2 * cosh_omega_2);
                real xc = 2.0 / (1.0 + exp_omega * exp_omega);

                // If quadrature node is sufficiently close to -1 or 1, declare it to be exactly
                // -1 or 1 and accumulate the quad weights
                volatile real test = real(1) - xc;
                if (abs(test - real(1)) <= real(0.0))
                {
                    if (snappedEndpoint)
                    {
                        out[count-3] += w;
                        out[count-1] += w;
                    }
                    else
                    {
                        out[count++] = -1.0;
                        out[count++] = w;
                        out[count++] = 1.0;
                        out[count++] = w;
                        snappedEndpoint = true;
                    }
                }
                else
                {
                    assert(!snappedEndpoint);
                    out[count++] = -1.0 + xc;
                    out[count++] = w;
                    out[count++] = 1.0 - xc;
                    out[count++] = w;
                }
            }

            assert(count % 2 == 0 && (snappedEndpoint && count <= 2*n || !snappedEndpoint && count == 2*n));

            // Normalise quadrature weights
            real sum = 0.0;
            for (int i = 0; i < count/2; ++i)
                sum += out[2*i + 1];
            real alpha = 2.0 / sum;
            for (int i = 0; i < count/2; ++i)
                out[2*i + 1] *= alpha;

            return count/2;
        }
    };
} // namespace algoim

#endif
