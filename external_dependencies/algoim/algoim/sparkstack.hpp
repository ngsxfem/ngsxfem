#ifndef ALGOIM_SPARKSTACK_HPP
#define ALGOIM_SPARKSTACK_HPP

// algoim::SparkStack<T> implements a fast, thread-safe, stack-based allocator,
// similar in function to alloca() but with additional guarantees regarding
// portability, alignment, and type consistency.

#include <vector>
#include "uvector.hpp"

namespace algoim
{
    template<typename T, int N>
    class xarray;

    template<typename T>
    class SparkStack
    {
        static constexpr size_t capacity = 1u << 23;
        static constexpr int capacity_line = __LINE__ - 1;

        template<typename ...R>
        static size_t alloc(T** ptr, size_t len, R... rest)
        {
            if (pos() + len > capacity)
            {
                std::cerr << "SparkStack<T = " << typeid(T).name() << ">: capacity=" << capacity << " and pos=" << pos() << " insufficient for request len=" << len << '\n';
                std::cerr << "    consider increasing const 'capacity', defined on line " << capacity_line << " in file " << __FILE__ << '\n';
                throw std::bad_alloc();
            }
            *ptr = base() + pos();
            pos() += len;
            if constexpr (sizeof...(rest) == 0)
                return len;
            else
                return len + alloc(rest...);
        }

        static T* base()
        {
            static thread_local std::vector<T> buff(capacity);
            return buff.data();
        }

        static ptrdiff_t& pos()
        {
            static thread_local ptrdiff_t pos_ = 0;
            return pos_;
        };

        size_t len_;

        SparkStack(const SparkStack&) = delete;
        SparkStack(SparkStack&&) = delete;
        SparkStack& operator=(const SparkStack&) = delete;
        SparkStack& operator=(SparkStack&&) = delete;

    public:

        // With parameters x0, n0, x1, n1, x2, n2, ..., allocate n0 elements and assign to x0, etc.
        template<typename ...R>
        explicit SparkStack(T** ptr, size_t len, R&&... rest)
        {
            len_ = alloc(ptr, len, rest...);
        }

        // With parameters value, x0, n0, x1, n1, x2, n2, ..., allocate n0 elements and assign to x0, ...,
        // and assign the given value to all n0*n1*n2*... values allocated
        template<typename ...R>
        explicit SparkStack(T value, T** ptr, size_t len, R&&... rest)
        {
            T* start = base() + pos();
            len_ = alloc(ptr, len, rest...);
            for (int i = 0; i < len_; ++i)
                *(start + i) = value;
        }

        // For each i, allocate ext(i) elements and assign to ptr(i)
        template<int N>
        explicit SparkStack(uvector<T*,N>& ptr, const uvector<int,N>& ext)
        {
            len_ = 0;
            for (int i = 0; i < N; ++i)
                len_ += alloc(&ptr(i), ext(i));
        }

        // Allocate enough elements for one or more xarray's having pre-set extent
        template<int ...N>
        explicit SparkStack(xarray<T,N>&... a)
        {
            len_ = (alloc(&a.data_, a.size()) + ...);
        }

        // Release memory when the SparkStack object goes out of scope
        ~SparkStack()
        {
            pos() -= len_;
        }
    };

    #define algoim_CONCAT2(x, y) x ## y
    #define algoim_CONCAT(x, y) algoim_CONCAT2(x, y)
    #define algoim_spark_alloc(T, ...) SparkStack<T> algoim_CONCAT(spark_alloc_var_, __LINE__)(__VA_ARGS__)
    #define algoim_spark_alloc_def(T, val, ...) SparkStack<T> algoim_CONCAT(spark_alloc_var_, __LINE__)(val, __VA_ARGS__)
    #define algoim_spark_alloc_vec(T, ptr, ext) SparkStack<T> algoim_CONCAT(spark_alloc_var_, __LINE__)(ptr, ext)

} // namespace algoim

#endif
