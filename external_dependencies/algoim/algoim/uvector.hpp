#ifndef ALGOIM_UVECTOR_HPP
#define ALGOIM_UVECTOR_HPP

// algoim::uvector implements a vector with compile-time known extent, including
// arithmetic expressions and various utility methods. It could be viewed as, e.g.,
// a light-weight, self-contained, re-implementation of blitz's TinyVector<T,N>,
// (https://en.wikipedia.org/wiki/Blitz%2B%2B), using more modern C++17/20 
// techniques.

#include <type_traits>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <cmath>

namespace algoim
{
    template<typename T, int N>
    class uvector;

    template<int N, typename T>
    struct uvector_expr
    {
        T eval;
    };

    namespace detail
    {
        template<typename T, int N>
        constexpr auto eval(const uvector<T,N>& e, int i)
        {
            return e(i);
        }

        template<int N, typename T>
        constexpr auto eval(const uvector_expr<N,T>& expr, int i)
        {
            return expr.eval(i);
        }

        template<typename E>
        constexpr auto eval(const E& e, int i)
        {
            return e;
        }

        template<int N, typename T>
        constexpr auto make_uvector_expr(const T& eval)
        {
            return uvector_expr<N,T>{eval};
        }

        template<typename T>
        struct extent_of { static constexpr int value = 0; };

        template<typename T, int N>
        struct extent_of<uvector<T,N>> { static constexpr int value = N; };

        template<int N, typename T>
        struct extent_of<uvector_expr<N,T>> { static constexpr int value = N; };

        template<typename T>
        using extent = std::integral_constant<int, extent_of<std::decay_t<T>>::value>;
    } // namespace detail

    // binary uvector expressions
    #define algoim_decl_binary_op(op)                                                                                           \
    template<typename E1, typename E2,                                                                                          \
        std::enable_if_t<(detail::extent<E1>::value > 0) || (detail::extent<E2>::value > 0), bool> = true>                      \
    auto operator op (E1&& e1, E2&& e2)                                                                                         \
    {                                                                                                                           \
        static_assert(detail::extent<E1>::value == 0 || detail::extent<E2>::value == 0 ||                                       \
            detail::extent<E1>::value == detail::extent<E2>::value,                                                             \
            "incompatible extents in algoim::uvector expression");                                                              \
        constexpr int N = std::max(detail::extent<E1>::value, detail::extent<E2>::value);                                       \
        if constexpr (std::is_lvalue_reference<E1&&>::value)                                                                    \
            if constexpr (std::is_lvalue_reference<E2&&>::value)                                                                \
                return detail::make_uvector_expr<N>([&e1,&e2](int i){ return detail::eval(e1, i) op detail::eval(e2, i); });    \
            else                                                                                                                \
                return detail::make_uvector_expr<N>([&e1, e2](int i){ return detail::eval(e1, i) op detail::eval(e2, i); });    \
        else                                                                                                                    \
            if constexpr (std::is_lvalue_reference<E2&&>::value)                                                                \
                return detail::make_uvector_expr<N>([ e1,&e2](int i){ return detail::eval(e1, i) op detail::eval(e2, i); });    \
            else                                                                                                                \
                return detail::make_uvector_expr<N>([ e1, e2](int i){ return detail::eval(e1, i) op detail::eval(e2, i); });    \
    }
    algoim_decl_binary_op(+)
    algoim_decl_binary_op(-)
    algoim_decl_binary_op(*)
    algoim_decl_binary_op(/)
    algoim_decl_binary_op(<)
    algoim_decl_binary_op(<=)
    algoim_decl_binary_op(>)
    algoim_decl_binary_op(>=)
    algoim_decl_binary_op(==)
    algoim_decl_binary_op(!=)
    algoim_decl_binary_op(%)
    #undef algoim_decl_binary_op

    // unary negation operator for uvector or uvector expressions
    template<typename E, std::enable_if_t<detail::extent<E>::value >= 1, bool> = true>
    auto operator- (const E& e)
    {
        constexpr int N = detail::extent<E>::value;
        if constexpr (std::is_lvalue_reference<E&&>::value)   
            return detail::make_uvector_expr<N>([&e](int i){ return -detail::eval(e, i); });
        else
            return detail::make_uvector_expr<N>([ e](int i){ return -detail::eval(e, i); });
    }

    // abs of uvector or uvector expressions
    template<typename E, std::enable_if_t<detail::extent<E>::value >= 1, bool> = true>
    auto abs (const E& e)
    {
        using std::abs;
        constexpr int N = detail::extent<E>::value;
        if constexpr (std::is_lvalue_reference<E&&>::value)   
            return detail::make_uvector_expr<N>([&e](int i){ return abs(detail::eval(e, i)); });
        else
            return detail::make_uvector_expr<N>([ e](int i){ return abs(detail::eval(e, i)); });
    }

    // uvector class
    template<typename T, int N>
    class uvector
    {
        T data_[N];
    public:
        // default constructor zero initialises
        uvector()
        {
            memset(data_, 0, sizeof(data_));
        }

        template<typename ...U, std::enable_if_t<N >= 2 && sizeof...(U) == N, bool> = true>
        explicit uvector(U&&... args) : data_{static_cast<T>(std::forward<U>(args))...}
        {
        }

        template<typename U>
        uvector(const U& x)
        {
            *this = x;
        }

              T& operator() (int i)       { return data_[i]; }
        const T& operator() (int i) const { return data_[i]; }

              T* data()       { return data_; }
        const T* data() const { return data_; }

        // assignment and update operators, e.g., operator=, operator+=, operator*=, etc.
        #define algoim_decl_assign_op(op)                                                   \
        template<typename E>                                                                \
        uvector& operator op (const E& x)                                                   \
        {                                                                                   \
            static_assert(detail::extent<E>::value == 0 || detail::extent<E>::value == N,   \
                "incompatible extents in algoim::uvector expression");                      \
            for (int i = 0; i < N; ++i)                                                     \
                data_[i] op detail::eval(x, i);                                             \
            return *this;                                                                   \
        }
        algoim_decl_assign_op(=)
        algoim_decl_assign_op(+=)
        algoim_decl_assign_op(-=)
        algoim_decl_assign_op(*=)
        algoim_decl_assign_op(/=)
        #undef algoim_decl_assign_op
    };

    // some methods, particularly those which remove components of a vector, benefit
    // from having a zero-length uvector; it is simply an empty class
    template<typename T>
    class uvector<T,0>
    {
    };

    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    std::ostream& operator<< (std::ostream& o, const T& u)
    {
        constexpr int N = detail::extent<T>::value;
        o << '(';
        for (int i = 0; i < N; ++i)
            o << detail::eval(u, i) << (i + 1 < N ? ',' : ')');
        return o;
    }

    // minimum component of a uvector or uvector expression
    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    auto min(const T& u)
    {
        constexpr int N = detail::extent<T>::value;
        auto x = detail::eval(u, 0);
        for (int i = 1; i < N; ++i)
        {
            auto y = detail::eval(u, i);
            if (y < x)
                x = y;
        }
        return x;
    }

    // maximum component of a uvector or uvector expression
    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    auto max(const T& u)
    {
        constexpr int N = detail::extent<T>::value;
        auto x = detail::eval(u, 0);
        for (int i = 1; i < N; ++i)
        {
            auto y = detail::eval(u, i);
            if (y > x)
                x = y;
        }
        return x;
    }

    // index where the minimum of a uvector or uvector expression is attained; if the minimum
    // is attained at multiple places, the smallest corresponding index is returned
    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    int argmin(const T& u)
    {
        constexpr int N = detail::extent<T>::value;
        auto x = detail::eval(u, 0);
        int ind = 0;
        for (int i = 1; i < N; ++i)
        {
            auto y = detail::eval(u, i);
            if (y < x)
            {
                x = y;
                ind = i;
            }
        }
        return ind;
    }

    // index where the maximum of a uvector or uvector expression is attained; if the maximum
    // is attained at multiple places, the smallest corresponding index is returned
    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    int argmax(const T& u)
    {
        constexpr int N = detail::extent<T>::value;
        auto x = detail::eval(u, 0);
        int ind = 0;
        for (int i = 1; i < N; ++i)
        {
            auto y = detail::eval(u, i);
            if (y > x)
            {
                x = y;
                ind = i;
            }
        }
        return ind;
    }

    // sum of components of a uvector or uvector expression
    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    auto sum(const T& u)
    {
        constexpr int N = detail::extent<T>::value;
        auto x = detail::eval(u, 0);
        for (int i = 1; i < N; ++i)
            x += detail::eval(u, i);
        return x;
    }

    // product of components of a uvector or uvector expression
    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    auto prod(const T& u)
    {
        constexpr int N = detail::extent<T>::value;
        auto x = detail::eval(u, 0);
        for (int i = 1; i < N; ++i)
            x *= detail::eval(u, i);
        return x;
    }

    // product of components of a uvector or uvector expression, excluding
    // a given index; the empty product returns one
    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    auto prod(const T& u, int excluded_index)
    {
        constexpr int N = detail::extent<T>::value;
        auto x = detail::eval(u, 0); // deduce return type
        x = 1;
        for (int i = 0; i < N; ++i)
            if (i != excluded_index)
                x *= detail::eval(u, i);
        return x;
    }

    // dot product of two uvector or uvector expressions
    template<typename E1, typename E2, std::enable_if_t<detail::extent<E1>::value == detail::extent<E2>::value && detail::extent<E1>::value >= 1, bool> = true>
    auto dot(const E1& u, const E2& v)
    {
        constexpr int N = detail::extent<E1>::value;
        auto x = detail::eval(u, 0) * detail::eval(v, 0);
        for (int i = 1; i < N; ++i)
            x += detail::eval(u, i) * detail::eval(v, i);
        return x;
    }

    // component-wise max of two uvectors
    template<typename T, int N>
    auto max(const uvector<T,N>& u, const uvector<T,N>& v)
    {
        using std::max;
        uvector<T,N> x;
        for (int i = 0; i < N; ++i)
            x(i) = max(u(i), v(i));
        return x;
    }

    // returns false iff all components of u yield false when typecast to a bool
    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    bool any(const T& u)
    {
        constexpr int N = detail::extent<T>::value;
        bool x = false;
        for (int i = 0; i < N; ++i)
            x |= static_cast<bool>(detail::eval(u, i));
        return x;
    }

    // returns true iff all components of u yield true when typecast to a bool
    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    bool all(const T& u)
    {
        constexpr int N = detail::extent<T>::value;
        bool x = true;
        for (int i = 0; i < N; ++i)
            x &= static_cast<bool>(detail::eval(u, i));
        return x;
    }

    // squared Euclidean norm of a uvector or uvector expression
    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    auto sqrnorm(const T& u)
    {
        return dot(u, u);
    }

    // Euclidean norm of a uvector or uvector expression
    template<typename T, std::enable_if_t<detail::extent<T>::value >= 1, bool> = true>
    auto norm(const T& u)
    {
        using std::sqrt;
        return sqrt(sqrnorm(u));
    }

    // alter the extent of a uvector by truncation or zero backfill
    template<int M, typename T, int N>
    uvector<T,M> alter_extent(const uvector<T,N>& u)
    {
        uvector<T,M> v;
        for (int i = 0; i < (M < N ? M : N); ++i) // truncate if M < N
            v(i) = u(i);
        for (int i = N; i < M; ++i) // backfill with zeros if M > N
            memset(&v(i), 0, sizeof(v(i)));
        return v;
    }

    // remove a component from a uvector
    template<typename T, int N>
    uvector<T,N-1> remove_component(const uvector<T,N>& u, int index)
    {
        uvector<T,N-1> v;
        if constexpr (N > 1)
            for (int i = 0; i < N - 1; ++i)
                v(i) = u(i < index ? i : i + 1);
        return v;
    }

    // add a component to a vector, consistent with remove_component semantics
    template<typename T, int N>
    uvector<T,N+1> add_component(const uvector<T,N>& u, int index, T value)
    {
        if constexpr (N == 0)
            return uvector<T,1>(value);
        else
        {
            uvector<T,N+1> v;
            for (int i = 0; i < N + 1; ++i)
                v(i) = i < index ? u(i) : (i == index ? value : u(i - 1));
            return v;
        }
    }

    // set the i'th component of a uvector to a given value
    template<typename T, int N>
    uvector<T,N> set_component(uvector<T,N> u, int i, T value)
    {
        u(i) = value;
        return u;
    }

    // increment the i'th component of a uvector by a given value
    template<typename T, int N>
    uvector<T,N> inc_component(uvector<T,N> u, int i, T value)
    {
        u(i) += value;
        return u;
    }
} // namespace algoim

#endif
