
#pragma once

#include <cmath>
#include <type_traits>
#include <ostream>

#include "common/math/Math.hpp"

#ifndef HD
    #ifdef __CUDACC__
        #define HD __host__ __device__
    #else
        #define HD
    #endif
#endif

template <int N, typename Real = float>
struct VecN 
{
    static_assert(N > 0, "VecN dimension must be positive.");
    static_assert(
        std::is_floating_point<Real>::value,
        "VecN<Real> requires a floating-point scalar type.");

    static constexpr int Size = N;

    Real v[N];

    template <int K>
    struct SegmentProxy
    {
        static_assert(K >= 1 && K <= N,
                    "Segment size must be between 1 and N");

        HD SegmentProxy(VecN& parent, int start)
            : m_parent(parent), m_start(start)
        {
        }

        // Read/write component access.
        HD Real& operator[](int i)
        {
            return m_parent[m_start + i];
        }

        HD const Real& operator[](int i) const
        {
            return m_parent[m_start + i];
        }

        // Convert the segment to an independent vector when needed.
        HD operator VecN<K, Real>() const
        {
            VecN<K, Real> result;
            for (int i = 0; i < K; ++i)
                result[i] = m_parent[m_start + i];

            return result;
        }

        // Assign from a regular vector.
        HD SegmentProxy& operator=(const VecN<K, Real>& rhs)
        {
            for (int i = 0; i < K; ++i)
                m_parent[m_start + i] = rhs[i];

            return *this;
        }

        // Assign from another segment proxy.
        HD SegmentProxy& operator=(const SegmentProxy& rhs)
        {
            // Copy through a temporary to handle overlapping segments safely.
            VecN<K, Real> tmp = static_cast<VecN<K, Real>>(rhs);
            return *this = tmp;
        }

    private:
        VecN& m_parent;
        int m_start;
    };

    // Constructors

    // Default construction leaves components uninitialized.
    HD VecN() = default;

    HD explicit VecN(Real s) {
        for (int i = 0; i < N; ++i)
            v[i] = s;
    }

    // Construct from exactly N scalar arguments.
    // Available when N == 2.
    template <int M = N, typename = typename
        std::enable_if<M == 2>::type>
    HD VecN(Real x, Real y) : v{x, y} {}

    // Construct from exactly N scalar arguments.
    // Available when N == 3.
    template <int M = N, typename = typename
        std::enable_if<M == 3>::type>
    HD VecN(Real x, Real y, Real z) : v{x, y, z} {}

    // Construct from exactly N scalar arguments.
    // Available when N == 4.
    template <int M = N, typename = typename
        std::enable_if<M == 4>::type>
    HD VecN(Real x, Real y, Real z, Real w) : v{x, y, z, w} {}

    // Construct from exactly N scalar arguments.
    // Available when N == 6.
    template <int M = N, typename = typename
        std::enable_if<M == 6>::type>
    HD VecN(Real x1, Real x2, Real x3, Real x4, Real x5, Real x6) 
        : v{x1, x2, x3, x4, x5, x6} {}

    // Construct by converting another vector type.
    template <typename OtherReal>
    HD explicit VecN(const VecN<N, OtherReal>& other) {
        for (int i = 0; i < N; ++i)
            v[i] = static_cast<Real>(other[i]);
    }

    // Indexed access (unchecked).
    HD Real& operator[](int i) {
        return v[i];
    }

    HD const Real& operator[](int i) const {
        return v[i];
    }

    HD Real& operator()(int i) {
        return v[i];
    }

    HD const Real& operator()(int i) const {
        return v[i];
    }

    HD const Real* data() const { return v; }
    HD Real* data() { return v; }

    // Named accessors for common dimensions.
    // These avoid assuming that array elements are separate
    // struct members, unlike pointer arithmetic from &v[0].
    HD Real& x() { return v[0]; }
    HD const Real& x() const { return v[0]; }

    HD Real& y() { return v[1]; }
    HD const Real& y() const { return v[1]; }

    HD Real& z() { return v[2]; }
    HD const Real& z() const { return v[2]; }

    HD Real& w() { return v[3]; }
    HD const Real& w() const { return v[3]; }

    // First K components: writable proxy.
    template <int K>
    HD SegmentProxy<K> head()
    {
        static_assert(K >= 1 && K <= N,
                    "head<K>: K must be between 1 and N");
        return SegmentProxy<K>(*this, 0);
    }

    // First K components: read-only copy.
    template <int K>
    HD VecN<K, Real> head() const
    {
        static_assert(K >= 1 && K <= N,
                    "head<K>: K must be between 1 and N");

        VecN<K, Real> result;
        for (int i = 0; i < K; ++i)
            result[i] = v[i];

        return result;
    }

    // Last K components: writable proxy.
    template <int K>
    HD SegmentProxy<K> tail()
    {
        static_assert(K >= 1 && K <= N,
                    "tail<K>: K must be between 1 and N");
        return SegmentProxy<K>(*this, N - K);
    }

    // Last K components: read-only copy.
    template <int K>
    HD VecN<K, Real> tail() const
    {
        static_assert(K >= 1 && K <= N,
                    "tail<K>: K must be between 1 and N");

        VecN<K, Real> result;
        for (int i = 0; i < K; ++i)
            result[i] = v[N - K + i];

        return result;
    }

    // Factories

    HD static VecN Zero() {
        return VecN(Real(0));
    }

    HD static VecN Ones() {
        return VecN(Real(1));
    }

    // Unary operators

    HD VecN operator+() const {
        return *this;
    }

    HD VecN operator-() const {
        VecN result;
        for (int i = 0; i < N; ++i)
            result[i] = -v[i];
        return result;
    }

    // Vector arithmetic

    HD VecN operator+(const VecN& b) const {
        VecN result;
        for (int i = 0; i < N; ++i)
            result[i] = v[i] + b[i];
        return result;
    }

    HD VecN operator-(const VecN& b) const {
        VecN result;
        for (int i = 0; i < N; ++i)
            result[i] = v[i] - b[i];
        return result;
    }

    HD VecN& operator+=(const VecN& b) {
        for (int i = 0; i < N; ++i)
            v[i] += b[i];
        return *this;
    }

    HD VecN& operator-=(const VecN& b) {
        for (int i = 0; i < N; ++i)
            v[i] -= b[i];
        return *this;
    }

    // Scalar arithmetic

    template <
        typename Scalar,
        typename = typename std::enable_if<
            std::is_arithmetic<Scalar>::value
        >::type
    >
    HD VecN operator*(Scalar s) const {
        VecN result;
        for (int i = 0; i < N; ++i)
            result[i] = v[i] * s;
        return result;
    }

    template <
        typename Scalar,
        typename = typename std::enable_if<
            std::is_arithmetic<Scalar>::value
        >::type
    >
    HD VecN operator/(Scalar s) const {
        const Real inv = Real(1) / s;
        return *this * inv;
    }

    template <
        typename Scalar,
        typename = typename std::enable_if<
            std::is_arithmetic<Scalar>::value
        >::type
    >
    HD VecN& operator*=(Scalar s) {
        for (int i = 0; i < N; ++i)
            v[i] *= s;
        return *this;
    }

    template <
        typename Scalar,
        typename = typename std::enable_if<
            std::is_arithmetic<Scalar>::value
        >::type
    >
    HD VecN& operator/=(Scalar s) {
        const Real inv = Real(1) / s;
        return *this *= inv;
    }

    template <
        typename Scalar,
        typename = typename std::enable_if<
            std::is_arithmetic<Scalar>::value
        >::type
    >
    HD VecN operator+(Scalar s) const {
        VecN result;
        for (int i = 0; i < N; ++i)
            result[i] = v[i] + s;
        return result;
    }

    template <
        typename Scalar,
        typename = typename std::enable_if<
            std::is_arithmetic<Scalar>::value
        >::type
    >
    HD VecN operator-(Scalar s) const {
        VecN result;
        for (int i = 0; i < N; ++i)
            result[i] = v[i] - s;
        return result;
    }

    template <
        typename Scalar,
        typename = typename std::enable_if<
            std::is_arithmetic<Scalar>::value
        >::type
    >
    HD VecN& operator+=(Scalar s) {
        for (int i = 0; i < N; ++i)
            v[i] += s;
        return *this;
    }

    template <
        typename Scalar,
        typename = typename std::enable_if<
            std::is_arithmetic<Scalar>::value
        >::type
    >
    HD VecN& operator-=(Scalar s) {
        for (int i = 0; i < N; ++i)
            v[i] -= s;
        return *this;
    }

    // Dot product

    HD Real dot(const VecN& b) const {
        Real result = Real(0);
        for (int i = 0; i < N; ++i)
            result += v[i] * b[i];
        return result;
    }

    // Cross product
    template <int M = N, typename = typename
        std::enable_if<M == 3>::type>
    HD VecN<3, Real> cross(const VecN<3, Real>& b) const
    {
        return VecN<3,Real>(
            y() * b.z() - z() * b.y(),
            z() * b.x() - x() * b.z(),
            x() * b.y() - y() * b.x()
        );
    }

    // Norms

    HD Real squaredNorm() const {
        return dot(*this);
    }

    HD Real norm() const {
        using ::sqrt;
        return static_cast<Real>(sqrt(squaredNorm()));
    }

    // Normalize in place.
    // Returns false if the norm is too small or non-finite.
    HD bool normalize(Real eps = Real(1e-8)) {
        const Real n2 = squaredNorm();

        if (!isfinite(n2) || n2 <= eps * eps)
            return false;

        using ::sqrt;
        const Real inv = Real(1) / static_cast<Real>(sqrt(n2));

        *this *= inv;
        return true;
    }

    // Normalized copy; returns zero if normalization fails.
    HD VecN normalized(Real eps = Real(1e-8)) const {
        VecN result = *this;
        if (!result.normalize(eps))
            return Zero();
        return result;
    }

    // Component-wise operations

    HD VecN cwiseProduct(const VecN& b) const {
        VecN result;
        for (int i = 0; i < N; ++i)
            result[i] = v[i] * b[i];
        return result;
    }

    HD VecN cwiseQuotient(const VecN& b) const {
        VecN result;
        for (int i = 0; i < N; ++i)
            result[i] = v[i] / b[i];
        return result;
    }

    HD VecN cwiseMin(const VecN& b) const {
        VecN result;
        for (int i = 0; i < N; ++i)
            result[i] = v[i] < b[i] ? v[i] : b[i];
        return result;
    }

    HD VecN cwiseMax(const VecN& b) const {
        VecN result;
        for (int i = 0; i < N; ++i)
            result[i] = v[i] > b[i] ? v[i] : b[i];
        return result;
    }

    HD VecN cwiseAbs() const {
        VecN result;
        for (int i = 0; i < N; ++i) {
            using ::fabs;
            result[i] = static_cast<Real>(fabs(v[i]));
        }
        return result;
    }

    HD Real minCoeff() const {
        Real result = v[0];
        for (int i = 1; i < N; ++i)
            result = v[i] < result ? v[i] : result;
        return result;
    }

    HD Real maxCoeff() const {
        Real result = v[0];
        for (int i = 1; i < N; ++i)
            result = v[i] > result ? v[i] : result;
        return result;
    }

    // Distance

    HD Real squaredDistance(const VecN& b) const {
        return (*this - b).squaredNorm();
    }

    HD Real distance(const VecN& b) const {
        using ::sqrt;
        return static_cast<Real>(sqrt(squaredDistance(b)));
    }

    // Exact equality; prefer a tolerance for floating-point comparisons.
    HD bool operator==(const VecN& b) const {
        for (int i = 0; i < N; ++i) {
            if (v[i] != b[i])
                return false;
        }
        return true;
    }

    HD bool operator!=(const VecN& b) const {
        return !(*this == b);
    }

    // Orthogonal vector
    HD VecN<N, Real> unitOrthogonal() const
    {
        VecN<N, Real> result = VecN<N, Real>::Zero();

        if constexpr (N == 2)
        {
            // Rotate 90 degrees.
            result[0] = -v[1];
            result[1] =  v[0];
        }
        else if constexpr (N == 3)
        {
            // Choose the coordinate axis corresponding to the
            // smallest absolute component to improve numerical stability.
            const Real ax = absValue(v[0]);
            const Real ay = absValue(v[1]);
            const Real az = absValue(v[2]);

            if (ax < ay)
            {
                if (ax < az)
                {
                    result[1] = -v[2];
                    result[2] =  v[1];
                }
                else
                {
                    result[0] = -v[1];
                    result[1] =  v[0];
                }
            }
            else
            {
                if (ay < az)
                {
                    result[0] = -v[2];
                    result[2] =  v[0];
                }
                else
                {
                    result[0] = -v[1];
                    result[1] =  v[0];
                }
            }
        }
        else
        {
            // General N-dimensional fallback:
            // choose the coordinate axis least aligned with this vector,
            // project that axis onto the orthogonal complement, then normalize.
            int k = 0;
            for (int i = 1; i < N; ++i)
            {
                if (absValue(v[i]) < absValue(v[k]))
                    k = i;
            }

            const Real n2 = squaredNorm();
            if (n2 <= Real(0))
                return result;

            result[k] = Real(1);
            const Real projection = v[k] / n2;

            for (int i = 0; i < N; ++i)
                result[i] -= projection * v[i];
        }

        result.normalize();
        return result;
    }
};

// Scalar multiplication with scalar on the left.
template <int N, typename Real, typename Scalar,
        typename = typename std::enable_if<
            std::is_arithmetic<Scalar>::value
        >::type
    >
HD inline VecN<N, Real> operator*(Scalar s, const VecN<N, Real>& a) {
    return a * s;
}

// output stream
template <int N, typename Real>
std::ostream& operator<<(std::ostream& os, const VecN<N, Real>& v)
{
    os << "[";
    for (int i = 0; i < N; ++i)
    {
        if (i > 0)
            os << ", ";
        os << v[i];
    }
    os << "]";
    return os;
}

// Common aliases.
// Vec3 and Vec6 retain the float-based names from the earlier classes.
template <typename Real = float>
using Vec2 = VecN<2, Real>;

template <typename Real = float>
using Vec3 = VecN<3, Real>;

template <typename Real = float>
using Vec4 = VecN<4, Real>;

template <typename Real = float>
using Vec6 = VecN<6, Real>;

