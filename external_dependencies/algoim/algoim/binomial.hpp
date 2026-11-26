#ifndef ALGOIM_BINOMIAL_HPP
#define ALGOIM_BINOMIAL_HPP

// algoim::Binomial implements a few simple methods to calculate binomial
// coefficients. For especially large coefficients, there are better methods
// (faster, more stable) than the ones used here, but these are sufficient
// for their application in the Algoim library

#include <array>
#include <vector>
#include <unordered_map>
#include "real.hpp"

namespace algoim
{
    struct Binomial
    {
        // Compute a row of binomial coefficients binom(n,i), i = 0,...,n;
        // it is assumed that 'out' has length at least n+1
        static void compute_row(int n, real* out)
        {
            out[0] = 1;
            if (n == 0)
                return;
            out[1] = n;
            for (int k = 2; k <= n/2; ++k)
                out[k] = (out[k - 1] * (n + 1 - k)) / k;
            for (int k = 0; k <= n/2; ++k)
                out[n-k] = out[k];
        }

        // Compute (and cache) all binomial coefficients binom(n,i), i = 0,...,n;
        // for n small enough, the row is obtained directly from a dense table;
        // for n large, the row is computed and cached in a sparse hash table
        static const real* row(int n)
        {
            static constexpr int m = 31;
            static const auto precomputed = []()
            {
                std::array<real,((m+1)*(m+2))/2> d;
                int i = 0;
                for (int n = 0; n <= m; ++n)
                {
                    compute_row(n, &d[i]);
                    i += n + 1;
                }
                return d;
            }();

            if (n <= m)
                return &precomputed[ (n*(n+1))/2 ];

            // All other n are stored in a thread_local hash table
            static thread_local std::unordered_map<int,std::vector<real>> coeff;
            auto& b = coeff[n];
            if (b.empty())
            {
                b.resize(n + 1);
                compute_row(n, b.data());
            }
            return b.data();
        };

        // returns the binomial coefficient (n \\ i)
        static real c(int n, int i)
        {
            return row(n)[i];
        }
    };
} // namespace algoim

#endif
