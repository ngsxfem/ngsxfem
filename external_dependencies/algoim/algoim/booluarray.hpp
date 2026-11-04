#ifndef ALGOIM_BOOLUARRAY_HPP
#define ALGOIM_BOOLUARRAY_HPP

// algoim::booluarray

#include <bitset>
#include "uvector.hpp"

namespace algoim
{
    namespace booluarray_detail
    {
        constexpr int pow(int base, int exp)
        {
            return exp == 0 ? 1 : base * pow(base, exp - 1);
        }

        template<int E, int N>
        constexpr int furl(const uvector<int,N>& i)
        {
            int ind = i(0);
            for (int j = 1; j < N; ++j)
                ind = ind * E + i(j);
            return ind;
        }
    }

    // booluarray implements a simple N-dimensional array of booleans, with
    // compile-time extent E across all dimensions; it is essentially a basic
    // specialisation of uarray<bool,N,E>
    template<int N, int E>
    class booluarray
    {
        constexpr static int size = booluarray_detail::pow(E, N);
        std::bitset<size> bits;
    public:
        booluarray() {}

        booluarray(bool val)
        {
            if (val)
                bits.set();
        }

        bool operator() (const uvector<int,N>& i) const
        {
            return bits[booluarray_detail::furl<E>(i)];
        }

        auto operator() (const uvector<int,N>& i)
        {
            return bits[booluarray_detail::furl<E>(i)];
        }

        // returns true iff the entire array is false
        bool none() const
        {
            return bits.none();
        }
    };
} // namespace algoim

#endif
