#ifndef ALGOIM_XARRAY_HPP
#define ALGOIM_XARRAY_HPP

// algoim::xarray implements an N-dimensional array view of a memory
// block controlled by the user

#include <cassert>
#include "uvector.hpp"
#include "sparkstack.hpp"
#include "utility.hpp"

namespace algoim
{
    // MiniLoop<N> is an optimised version of MultiLoop<N> with zero lowerbound
    template<int N>
    class MiniLoop
    {
        uvector<int,N> i;
        int iexp;
        const uvector<int,N> ext;
    public:
        explicit MiniLoop(const uvector<int,N>& ext) : ext(ext), i(0), iexp(0) {}

        MiniLoop& operator++()
        {
            ++iexp;
            for (int dim = N - 1; dim >= 0; --dim)
            {
                if (++i(dim) < ext(dim))
                    return *this;
                if (dim == 0)
                    return *this;
                i(dim) = 0;
            }
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
            return i(0) < ext(0);
        }

        int furl() const
        {
            return iexp;
        }

        uvector<int,N> shifted(int dim, int amount) const
        {
            uvector<int,N> j = i;
            j(dim) += amount;
            return j;
        }
    };

    template<typename T>
    struct xarraySlice
    {
        T* ptr;
        int len;

        xarraySlice(xarraySlice&) = delete;
        xarraySlice(xarraySlice&&) = delete;

        xarraySlice& operator= (const T& x) { for (int i = 0; i < len; ++i) ptr[i] =  x; return *this; }
        xarraySlice& operator+=(const T& x) { for (int i = 0; i < len; ++i) ptr[i] += x; return *this; }
        xarraySlice& operator-=(const T& x) { for (int i = 0; i < len; ++i) ptr[i] -= x; return *this; }
        xarraySlice& operator*=(const T& x) { for (int i = 0; i < len; ++i) ptr[i] *= x; return *this; }
        xarraySlice& operator/=(const T& x) { for (int i = 0; i < len; ++i) ptr[i] /= x; return *this; }

        xarraySlice& operator= (const xarraySlice& x) { for (int i = 0; i < len; ++i) ptr[i] =  x.ptr[i]; return *this; }
        xarraySlice& operator+=(const xarraySlice& x) { for (int i = 0; i < len; ++i) ptr[i] += x.ptr[i]; return *this; }
        xarraySlice& operator-=(const xarraySlice& x) { for (int i = 0; i < len; ++i) ptr[i] -= x.ptr[i]; return *this; }
        xarraySlice& operator*=(const xarraySlice& x) { for (int i = 0; i < len; ++i) ptr[i] *= x.ptr[i]; return *this; }
        xarraySlice& operator/=(const xarraySlice& x) { for (int i = 0; i < len; ++i) ptr[i] /= x.ptr[i]; return *this; }

        template<typename X, typename Y>
        struct prod
        {
            const X& x;
            const Y& y;
        };

        xarraySlice& operator= (const prod<xarraySlice,T>& op) { for (int i = 0; i < len; ++i) ptr[i] =  op.x.ptr[i] * op.y; return *this; }
        xarraySlice& operator+=(const prod<xarraySlice,T>& op) { for (int i = 0; i < len; ++i) ptr[i] += op.x.ptr[i] * op.y; return *this; }
        xarraySlice& operator-=(const prod<xarraySlice,T>& op) { for (int i = 0; i < len; ++i) ptr[i] -= op.x.ptr[i] * op.y; return *this; }
        xarraySlice& operator*=(const prod<xarraySlice,T>& op) { for (int i = 0; i < len; ++i) ptr[i] *= op.x.ptr[i] * op.y; return *this; }
        xarraySlice& operator/=(const prod<xarraySlice,T>& op) { for (int i = 0; i < len; ++i) ptr[i] /= op.x.ptr[i] * op.y; return *this; }
    };

    template<typename T>
    auto operator* (const xarraySlice<T>& x, const T& y)
    {
        return typename xarraySlice<T>::template prod<xarraySlice<T>,T>{x, y};
    };

    template<typename T>
    void swap(const xarraySlice<T>& x, const xarraySlice<T>& y)
    {
        using std::swap;
        for (int i = 0; i < x.len; ++i)
            swap(x.ptr[i], y.ptr[i]);
    }

    // algoim::xarray implements an N-dimensional array view of a memory
    // block controlled by the user
    template<typename T, int N>
    class xarray
    {
        T* data_;
        uvector<int,N> ext_;
        friend class SparkStack<T>;
    public:
        // interpret the given block of memory as an N-dimensional array of the given extent
        xarray(T* data, const uvector<int,N>& ext) : data_(data), ext_(ext) {}

        xarray(const xarray&) = delete;

        xarray& operator=(const xarray& x)
        {
            assert(same_shape(x));
            for (int i = 0; i < size(); ++i)
                data_[i] = x.data_[i];
            return *this;
        }

        template<typename S>
        xarray& operator=(const S& x)
        {
            for (int i = 0; i < size(); ++i)
                data_[i] = x;
            return *this;
        }

        template<typename S>
        xarray& operator+=(const S& x)
        {
            for (int i = 0; i < size(); ++i)
                data_[i] += x;
            return *this;
        }

        xarray& operator+=(const xarray& x)
        {
            assert(same_shape(x));
            for (int i = 0; i < size(); ++i)
                data_[i] += x.data_[i];
            return *this;
        }

        template<typename S>
        xarray& operator-=(const S& x)
        {
            for (int i = 0; i < size(); ++i)
                data_[i] -= x;
            return *this;
        }

        xarray& operator-=(const xarray& x)
        {
            assert(same_shape(x));
            for (int i = 0; i < size(); ++i)
                data_[i] -= x.data_[i];
            return *this;
        }

        template<typename S>
        xarray& operator*=(const S& x)
        {
            for (int i = 0; i < size(); ++i)
                data_[i] *= x;
            return *this;
        }

        const T* data() const { return data_; }
        T* data() { return data_; }

        // Accessors for user-expanded index, 0 <= i < prod(ext)
        const T& operator[](int i) const { return data_[i]; }
        T& operator[](int i) { return data_[i]; }

        // Slice and flatten this(i,:)
        template<int NN = N, std::enable_if_t<NN == 1, bool> = true>
        T& a(int i)
        {
            return data_[i];
        }

        // Slice and flatten this(i,:)
        template<int NN = N, std::enable_if_t<NN == 1, bool> = true>
        const T& a(int i) const
        {
            return data_[i];
        }

        // Slice and flatten this(i,:)
        template<int NN = N, std::enable_if_t<(NN > 1), bool> = true>
        auto a(int i) const
        {
            int span = prod(ext_, 0);
            return xarraySlice<T>{data_ + i * span, span};
        }

        const T& operator()(int i, int j) const
        {
            static_assert(N == 2, "N == 2 required for integer pair access");
            return data_[i*ext(1) + j];
        }

        T& operator()(int i, int j)
        {
            static_assert(N == 2, "N == 2 required for integer pair access");
            return data_[i*ext(1) + j];
        }

        // Accessors by multi-index
        const T& m(const uvector<int,N>& i) const
        {
            return data_[util::furl(i, ext_)];
        }

        // Accessors by multi-index
        T& m(const uvector<int,N>& i)
        {
            return data_[util::furl(i, ext_)];
        }

        // Accessors by mini-loop, where it is assumed the MiniLoop object's
        // associated extent is identical to this xarray
        const T& l(const MiniLoop<N>& i) const
        {
            return data_[i.furl()];
        }

        T& l(const MiniLoop<N>& i)
        {
            return data_[i.furl()];
        }

        const uvector<int,N>& ext() const
        {
            return ext_;
        }

        int ext(int i) const
        {
            return ext_(i);
        }

        int size() const
        {
            return prod(ext_);
        }

        MiniLoop<N> loop() const
        {
            return MiniLoop<N>(ext_);
        }

        // xarray is strided so that the inner-most dimension (dim 0) is the slowest varying;
        // as such, one may view an N-dimensionl xarray as a length ext(0) vector of (N-1)-
        // dimensional xarrays. Two views of this(i,:) are provided:
        //   - flatten() returns a 2D view, of dimensions ext(0) by ext(1) * ... * ext(N - 1)
        //   - slice(i) returns an (N-1)-D array of extent (ext(1), ..., ext(N-1))
        auto flatten() const
        {
            return xarray<T,2>(data_, uvector<int,2>{ext_(0), prod(ext_, 0)});
        }

        // xarray is strided so that the inner-most dimension (dim 0) is the slowest varying;
        // as such, one may view an N-dimensionl xarray as a length ext(0) vector of (N-1)-
        // dimensional xarrays. Two views of this(i,:) are provided:
        //   - flatten() returns a 2D view, of dimensions ext(0) by ext(1) * ... * ext(N - 1)
        //   - slice(i) returns an (N-1)-D array of extent (ext(1), ..., ext(N-1))
        auto slice(int i) const
        {
            return xarray<T,N-1>(data_ + i * prod(ext_, 0), remove_component(ext_, 0));
        }

        bool same_shape(const xarray& x) const
        {
            return all(x.ext_ == ext_);
        }

        T maxNorm() const
        {
            assert(size() > 0);
            using std::abs;
            using std::max;
            T m = abs(data_[0]);
            for (int i = 1; i < size(); ++i)
                m = max(m, abs(data_[i]));
            return m;
        }

        T min() const
        {
            assert(size() > 0);
            using std::min;
            T m = data_[0];
            for (int i = 1; i < size(); ++i)
                m = min(m, data_[i]);
            return m;
        }

        T max() const
        {
            assert(size() > 0);
            using std::max;
            T m = data_[0];
            for (int i = 1; i < size(); ++i)
                m = max(m, data_[i]);
            return m;
        }

        // Alter the extent of the view of the existing memory block; usually this
        // only makes sense when the new extent is smaller than the existing extent,
        // i.e., prod(new extent) <= prod(old extent)
        void alterExtent(const uvector<int,N>& ext)
        {
            ext_ = ext;
        }

        // flatten() and slice() return temporary objects which cannot be bound to non-const references;
        // applying .ref() to the temporary object allows this to be done and is an assertion by the
        // user that doing so is safe (relating to the output of flatten() or slice() being destroyed
        // at the end of the full-expression)
        xarray<T,N>& ref()
        {
            return *this;
        }
    };
} // namespace algoim

#endif
