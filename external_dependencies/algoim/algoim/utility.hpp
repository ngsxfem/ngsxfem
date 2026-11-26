#ifndef ALGOIM_UTILITY_HPP
#define ALGOIM_UTILITY_HPP

// Minor utility methods used throughout Algoim

namespace algoim::util 
{
    static_assert(std::is_same_v<real, double>, "Warning: pi constant may require redefining when real != double");
    static constexpr double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816;

    // square of u
    template<typename T>
    constexpr auto sqr(T u)
    {
        return u*u;
    }

    // cube of u
    template<typename T>
    constexpr auto cube(T u)
    {
        return u*u*u;
    }

    // sign of u:
    //   -1, if u < 0,
    //   +1, if u > 0,
    //    0, otherwise
    template <typename T>
    constexpr int sign(T u) noexcept
    {
        return (T(0) < u) - (u < T(0));
    }

    // collapse an N-dimensional multi-index into a scalar integer according to its location
    // in an N-dimensional grid of the given extent, such that the lowest dimension iterates
    // the slowest, and the highest dimension corresponds to the inner-most loop, consistent
    // with MultiLoop semantics. For example,
    //    furl( {i,j}, {m,n} ) = i*n + j
    //    furl( {i,j,k}, {m,n,o} ) = i*n*o + j*o + k
    template<int N>
    int furl(const uvector<int,N>& i, const uvector<int,N>& ext)
    {
        int ind = i(0);
        for (int j = 1; j < N; ++j)
            ind = ext(j) * ind + i(j);
        return ind;
    }

    // compute a Givens rotation
    template<typename T>
    void givens_get(const T& a, const T& b, T& c, T& s)
    {
        using std::abs;
        using std::sqrt;
        if (b == 0.0)
        {
            c = 1.0;
            s = 0.0;
        }
        else if (abs(b) > abs(a))
        {
            T tmp = a / b;
            s = T(1) / sqrt(1.0 + tmp*tmp);
            c = tmp * s;
        }
        else
        {
            T tmp = b / a;
            c = T(1) / sqrt(1.0 + tmp*tmp);
            s = tmp * c;
        }
    };

    // apply a Givens rotation
    template<typename T>
    void givens_rotate(T& x, T& y, T c, T s)
    {
        T a = x, b = y;
        x =  c * a + s * b;
        y = -s * a + c * b;
    };
} // namespace algoim::util

#endif
