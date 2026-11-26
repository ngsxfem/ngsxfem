#ifndef ALGOIM_MULTILOOP_HPP
#define ALGOIM_MULTILOOP_HPP

// algoim::MultiLoop for writing N-dimensional nested for loops

#include "uvector.hpp"

namespace algoim
{
    // MultiLoop is essentially an N-dimensional iterator for looping over the
    // coordinates of a Cartesian grid having indices min(0) <= i < max(0),
    // min(1) <= j < max(1), min(2) <= k < max(2), etc. The ordering is such that
    // dimension N-1 is inner-most, i.e., iterates the fastest, while dimension
    // 0 is outer-most and iterates the slowest.
    template<int N>
    class MultiLoop
    {
        uvector<int,N> i;
        const uvector<int,N> min, max;
        bool valid;
    public:
        MultiLoop(const uvector<int,N>& min, const uvector<int,N>& max)
            : i(min), min(min), max(max), valid(all(min < max))
        {}

        MultiLoop& operator++()
        {
            for (int dim = N - 1; dim >= 0; --dim)
            {
                if (++i(dim) < max(dim))
                    return *this;
                i(dim) = min(dim);
            }
            valid = false;
            return *this;
        }

        const uvector<int,N>& operator()() const
        {
            return i;
        }

        int operator()(int index) const
        {
            return i(index);
        }

        bool operator~() const
        {
            return valid;
        }
    };
} // namespace algoim

#endif
