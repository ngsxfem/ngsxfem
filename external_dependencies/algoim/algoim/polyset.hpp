#ifndef ALGOIM_POLYSET_HPP
#define ALGOIM_POLYSET_HPP

// algoim::PolySet

#include <cassert>
#include <vector>
#include "booluarray.hpp"

namespace algoim
{
    // PolySet implements a simple container to hold one or more Bernstein polynomials
    // and their associated masks
    template<int N, int E>
    struct PolySet
    {
        struct Poly
        {
            uvector<int,N> ext;             // Degree/extent of polynomial
            size_t offset;                  // Offset into buffer, storing the xarray<Real,N> polynomial data
            booluarray<N,E> mask;           // Mask
        };
        std::vector<real> buff;             // Memory buffer containing polynomial data
        std::vector<Poly> items;            // Record of contained polynomials

        // Access polynomial by index
        xarray<real,N> poly(size_t ind)
        {
            assert(0 <= ind && ind < items.size());
            return xarray<real,N>(&buff[items[ind].offset], items[ind].ext);
        }

        // Access mask by index
        booluarray<N,E>& mask(size_t ind)
        {
            assert(0 <= ind && ind < items.size());
            return items[ind].mask;
        }

        // Add a polynomial/mask pair to the container
        void push_back(const xarray<real,N>& p, const booluarray<N,E>& m)
        {
            items.push_back({p.ext(), buff.size(), m});
            buff.resize(buff.size() + p.size());
            poly(items.size() - 1) = p;
        }

        size_t count() const
        {
            return items.size();
        }
    };
} // namespace algoim

#endif
