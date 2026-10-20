
#pragma once

#include "common/math/VecN.hpp"

#include <cmath>
#include <type_traits>
#include <ostream>

#ifndef HD
    #ifdef __CUDACC__
        #define HD __host__ __device__
    #else
        #define HD
    #endif
#endif

/////////////////////////////////////////////////////////////////////

template <int Rows, int Cols, typename Real = float>
struct MatN {
    static_assert(Rows > 0 && Cols > 0,
                  "Matrix dimensions must be positive.");
    static_assert(std::is_floating_point<Real>::value,
                  "MatN requires a floating-point scalar.");

    using Scalar = Real;
    static constexpr int NumRows = Rows;
    static constexpr int NumCols = Cols;
    static constexpr int NumElements = Rows * Cols;

    // Row-major storage.
    Real m[NumElements];


    // for stream-based (Eigen style) initialization of matrices
    class MatNInitializer
    {
    public:
        HD explicit MatNInitializer(MatN<Rows, Cols, Real>& matrix)
            : m_matrix(matrix), m_index(0)
        {
        }

        HD MatNInitializer& operator,(Real value)
        {
            assert(m_index < Rows * Cols);
            m_matrix(m_index / Cols, m_index % Cols) = value;
            ++m_index;
            return *this;
        }

        HD void initialize(Real value)
        {
            m_matrix(0, 0) = value;
            m_index = 1;
        }

    private:
        MatN<Rows, Cols, Real>& m_matrix;
        int m_index;
    };

    template <int BlockRows, int BlockCols>
    class MatNBlock
    {
    public:
        HD MatNBlock(MatN<Rows, Cols, Real>& parent,
                    int startRow, int startCol)
            : m_parent(parent), m_startRow(startRow), m_startCol(startCol)
        {
        }

        // Assign a block from a matrix-like value.
        HD MatNBlock& operator=(const MatN<BlockRows, BlockCols, Real>& rhs)
        {
            for (int i = 0; i < BlockRows; ++i)
                for (int j = 0; j < BlockCols; ++j)
                    m_parent(m_startRow + i, m_startCol + j) = rhs(i, j);

            return *this;
        }

        // Assign a vector to a block when it has one column or one row.
        template <
            int R = BlockRows,
            int C = BlockCols,
            typename = typename std::enable_if<
                (R == 1 || C == 1)
            >::type
        >
        HD MatNBlock& operator=(const VecN<(BlockRows == 1 ? BlockCols : BlockRows), Real>& rhs)
        {
            if constexpr (BlockCols == 1)
            {
                for (int i = 0; i < BlockRows; ++i)
                    m_parent(m_startRow + i, m_startCol) = rhs[i];
            }
            else // BlockRows == 1
            {
                for (int i = 0; i < BlockCols; ++i)
                    m_parent(m_startRow, m_startCol + i) = rhs[i];
            }

            return *this;
        }

        // Element access through the block.
        HD Real& operator()(int i, int j)
        {
            return m_parent(m_startRow + i, m_startCol + j);
        }

        HD const Real& operator()(int i, int j) const
        {
            return m_parent(m_startRow + i, m_startCol + j);
        }

    private:
        MatN<Rows, Cols, Real>& m_parent;
        int m_startRow;
        int m_startCol;
    };


    // Default construction leaves elements uninitialized.
    HD MatN() = default;

    // Diagonal constructor for square matrices.
    template <int R = Rows, int C = Cols,
              typename = typename std::enable_if<R == C>::type>
    HD explicit MatN(Real diagonal) {
        for (int r = 0; r < Rows; ++r)
            for (int c = 0; c < Cols; ++c)
                (*this)(r, c) = (r == c) ? diagonal : Real(0);
    }

    HD static MatN Zero() {
        MatN result;
        for (int i = 0; i < NumElements; ++i)
            result.m[i] = Real(0);
        return result;
    }

    template <int R = Rows, int C = Cols,
              typename = typename std::enable_if<R == C>::type>
    HD static MatN Identity() {
        return MatN(Real(1));
    }

    HD MatNInitializer operator<<(Real value)
    {
        MatNInitializer initializer(*this);
        initializer.initialize(value);
        return initializer;
    }

    // Unchecked element access.
    HD Real& operator()(int row, int col) {
        return m[row * Cols + col];
    }

    HD const Real& operator()(int row, int col) const {
        return m[row * Cols + col];
    }

    HD Real& operator[](int i) {
        return m[i];
    }

    HD const Real& operator[](int i) const {
        return m[i];
    }

    HD const Real* data() const { return m; }
    HD Real* data() { return m; }

    // Row and column access.
    HD VecN<Cols, Real> row(int r) const {
        VecN<Cols, Real> result;
        for (int c = 0; c < Cols; ++c)
            result[c] = (*this)(r, c);
        return result;
    }

    HD VecN<Rows, Real> col(int c) const {
        VecN<Rows, Real> result;
        for (int r = 0; r < Rows; ++r)
            result[r] = (*this)(r, c);
        return result;
    }

    HD void setRow(int r, const VecN<Cols, Real>& v) {
        for (int c = 0; c < Cols; ++c)
            (*this)(r, c) = v[c];
    }

    HD void setCol(int c, const VecN<Rows, Real>& v) {
        for (int r = 0; r < Rows; ++r)
            (*this)(r, c) = v[r];
    }

    // block access
    template <int BlockRows, int BlockCols>
    HD MatNBlock<BlockRows, BlockCols>
    block(int startRow, int startCol)
    {
        return MatNBlock<BlockRows, BlockCols>(
            *this, startRow, startCol);
    }

    // Unary operators.
    HD MatN operator-() const {
        MatN result;
        for (int i = 0; i < NumElements; ++i)
            result[i] = -m[i];
        return result;
    }

    // Matrix addition and subtraction.
    HD MatN operator+(const MatN& b) const {
        MatN result;
        for (int i = 0; i < NumElements; ++i)
            result[i] = m[i] + b[i];
        return result;
    }

    HD MatN operator-(const MatN& b) const {
        MatN result;
        for (int i = 0; i < NumElements; ++i)
            result[i] = m[i] - b[i];
        return result;
    }

    HD MatN& operator+=(const MatN& b) {
        for (int i = 0; i < NumElements; ++i)
            m[i] += b[i];
        return *this;
    }

    HD MatN& operator-=(const MatN& b) {
        for (int i = 0; i < NumElements; ++i)
            m[i] -= b[i];
        return *this;
    }

    // Scalar operations.
    template <typename Scalar>
    HD MatN operator*(Scalar s) const {
        MatN result;
        for (int i = 0; i < NumElements; ++i)
            result[i] = m[i] * s;
        return result;
    }

    template <typename Scalar>
    HD MatN operator/(Scalar s) const {
        return *this * (Real(1) / s);
    }

    template <typename Scalar>
    HD MatN& operator*=(Scalar s) {
        for (int i = 0; i < NumElements; ++i)
            m[i] *= s;
        return *this;
    }

    template <typename Scalar>
    HD MatN& operator/=(Scalar s) {
        return *this *= (Real(1) / s);
    }

    // Matrix-vector multiplication: A * columnVector.
    HD VecN<Rows, Real> operator*(
        const VecN<Cols, Real>& v) const
    {
        VecN<Rows, Real> result = VecN<Rows, Real>::Zero();

        for (int r = 0; r < Rows; ++r)
            for (int c = 0; c < Cols; ++c)
                result[r] += (*this)(r, c) * v[c];

        return result;
    }

    // Matrix multiplication:
    // MatN<R,C> * MatN<C,K> -> MatN<R,K>.
    template <int OtherCols>
    HD MatN<Rows, OtherCols, Real>
    operator*(const MatN<Cols, OtherCols, Real>& b) const
    {
        MatN<Rows, OtherCols, Real> result =
            MatN<Rows, OtherCols, Real>::Zero();

        for (int i = 0; i < Rows; ++i)
        {
            for (int k = 0; k < Cols; ++k)
            {
                const Real aik = (*this)(i, k);

                for (int j = 0; j < OtherCols; ++j)
                {
                    result(i, j) += aik * b(k, j);
                }
            }
        }

        return result;
    }

    // scales rows by elements of a vector
    HD MatN scaleRows(const VecN<Rows, Real>& v) const
    {
        MatN result;

        for (int i = 0; i < Rows; ++i)
            for (int j = 0; j < Cols; ++j)
                result(i, j) = (*this)(i, j) * v[i];

        return result;
    }

    // scales cols by elements of a vector
    HD MatN scaleCols(const VecN<Cols, Real>& v) const
    {
        MatN result;

        for (int i = 0; i < Rows; ++i)
            for (int j = 0; j < Cols; ++j)
                result(i, j) = (*this)(i, j) * v[j];

        return result;
    }

    // Transpose swaps row and column dimensions.
    HD MatN<Cols, Rows, Real> transpose() const {
        MatN<Cols, Rows, Real> result;

        for (int r = 0; r < Rows; ++r)
            for (int c = 0; c < Cols; ++c)
                result(c, r) = (*this)(r, c);

        return result;
    }

    // coefficient-wise absolute value
    HD MatN<Cols, Rows, Real> cwiseAbs() const {
        MatN<Cols, Rows, Real> result;
        for (int i = 0; i < NumElements; ++i) {
            using ::fabs;
            result[i] = static_cast<Real>(fabs(m[i]));
        }
        return result;
    }

    // Frobenius norm and squared norm.
    HD Real squaredNorm() const {
        Real result = Real(0);
        for (int i = 0; i < NumElements; ++i)
            result += m[i] * m[i];
        return result;
    }

    HD Real norm() const {
        using ::sqrt;
        return static_cast<Real>(sqrt(squaredNorm()));
    }

    // Square-matrix operations.
    template <int R = Rows, int C = Cols,
              typename = typename std::enable_if<R == C>::type>
    HD Real trace() const {
        Real result = Real(0);
        for (int i = 0; i < Rows; ++i)
            result += (*this)(i, i);
        return result;
    }

    template <int R = Rows, int C = Cols,
              typename = typename std::enable_if<R == C>::type>
    HD Real determinant(Real relativeEps = Real(1e-6)) const {
        Real a[Rows][Cols];
        Real scale = Real(0);
        int sign = 1;

        for (int r = 0; r < Rows; ++r) {
            for (int c = 0; c < Cols; ++c) {
                const Real value = (*this)(r, c);
                if (!isfinite(value))
                    return Real(0);

                a[r][c] = value;
                const Real av = fabs(value);
                if (av > scale)
                    scale = av;
            }
        }

        if (scale == Real(0))
            return Real(0);

        const Real threshold = relativeEps * scale;
        Real det = Real(1);

        for (int k = 0; k < Rows; ++k) {
            int pivotRow = k;
            Real pivotAbs = fabs(a[k][k]);

            for (int r = k + 1; r < Rows; ++r) {
                const Real candidate = fabs(a[r][k]);
                if (candidate > pivotAbs) {
                    pivotAbs = candidate;
                    pivotRow = r;
                }
            }

            if (pivotAbs <= threshold)
                return Real(0);

            if (pivotRow != k) {
                for (int c = 0; c < Cols; ++c) {
                    const Real tmp = a[k][c];
                    a[k][c] = a[pivotRow][c];
                    a[pivotRow][c] = tmp;
                }
                sign = -sign;
            }

            const Real pivot = a[k][k];
            det *= pivot;

            for (int r = k + 1; r < Rows; ++r) {
                const Real factor = a[r][k] / pivot;
                for (int c = k + 1; c < Cols; ++c)
                    a[r][c] -= factor * a[k][c];
            }
        }

        return sign < 0 ? -det : det;
    }

    // Solve A*x = b using Gaussian elimination with partial pivoting.
    // Returns false on singular/non-finite pivots.
    // x remains unchanged on failure.
    template <int R = Rows, int C = Cols,
              typename = typename std::enable_if<R == C>::type>
    HD bool solve(
        const VecN<Rows, Real>& b,
        VecN<Rows, Real>& x,
        Real relativeEps = Real(1e-6)) const
    {
        Real a[Rows][Cols];
        Real rhs[Rows];
        Real rowScale[Rows];

        for (int r = 0; r < Rows; ++r)
        {
            rowScale[r] = Real(0);
            rhs[r] = b[r];

            for (int c = 0; c < Cols; ++c)
            {
                const Real value = (*this)(r, c);
                if (!isfinite(value))
                    return false;

                a[r][c] = value;

                const Real av = fabs(value);
                if (av > rowScale[r])
                    rowScale[r] = av;
            }

            if (rowScale[r] == Real(0))
                return false;
        }

        for (int k = 0; k < Rows; ++k) {
            int pivotRow = k;
            Real bestRatio = Real(-1);

            for (int r = k; r < Rows; ++r)
            {
                const Real ratio =
                    fabs(a[r][k]) / rowScale[r];

                if (ratio > bestRatio)
                {
                    bestRatio = ratio;
                    pivotRow = r;
                }
            }

            const Real pivotAbs = fabs(a[pivotRow][k]);

            if (!isfinite(pivotAbs) ||
                pivotAbs <= relativeEps * rowScale[pivotRow])
            {
                return false;
            }

            if (pivotRow != k) {
                for (int c = k; c < Cols; ++c) {
                    const Real tmp = a[k][c];
                    a[k][c] = a[pivotRow][c];
                    a[pivotRow][c] = tmp;
                }

                const Real tmp = rhs[k];
                rhs[k] = rhs[pivotRow];
                rhs[pivotRow] = tmp;
                const Real tmpScale = rowScale[k];
                rowScale[k] = rowScale[pivotRow];
                rowScale[pivotRow] = tmpScale;
            }

            const Real pivot = a[k][k];

            for (int r = k + 1; r < Rows; ++r) {
                const Real factor = a[r][k] / pivot;
                a[r][k] = Real(0);

                for (int c = k + 1; c < Cols; ++c)
                    a[r][c] -= factor * a[k][c];

                rhs[r] -= factor * rhs[k];
            }
        }

        VecN<Rows, Real> result;

        for (int r = Rows - 1; r >= 0; --r) {
            Real value = rhs[r];

            for (int c = r + 1; c < Rows; ++c)
                value -= a[r][c] * result[c];

            const Real pivot = a[r][r];
            if (!isfinite(pivot) || fabs(pivot) <= relativeEps * rowScale[r])
            {
                return false;
            }

            result[r] = value / pivot;
            if (!isfinite(result[r]))
                return false;
        }

        x = result;
        return true;
    }

    // Inverse using repeated solves, one column at a time.
    template <int R = Rows, int C = Cols,
              typename = typename std::enable_if<R == C>::type>
    HD bool tryInverse(
        MatN& result,
        Real relativeEps = Real(1e-6)) const
    {
        MatN inverse = Zero();

        for (int c = 0; c < Cols; ++c) {
            VecN<Rows, Real> b = VecN<Rows, Real>::Zero();
            b[c] = Real(1);

            VecN<Rows, Real> x;
            if (!solve(b, x, relativeEps))
                return false;

            for (int r = 0; r < Rows; ++r)
                inverse(r, c) = x[r];
        }

        result = inverse;
        return true;
    }
};

template <int R, int C, typename Real, typename Scalar>
HD inline MatN<R, C, Real> operator*(
    Scalar s, const MatN<R, C, Real>& a)
{
    return a * s;
}

template <int Rows, int Cols, typename Real>
std::ostream& operator<<(std::ostream& os, const MatN<Rows, Cols, Real>& A)
{
    for (int i = 0; i < Rows; ++i)
    {
        os << "[";
        for (int j = 0; j < Cols; ++j)
        {
            if (j > 0)
                os << ", ";
            os << A(i, j);
        }
        os << "]";

        if (i + 1 < Rows)
            os << '\n';
    }

    return os;
}


/////////////////////////////////////////////////////////////////

template <int M, int N, typename Real>
HD MatN<M, N, Real>
outerProduct(const VecN<M, Real>& a,
             const VecN<N, Real>& b)
{
    MatN<M, N, Real> result = MatN<M, N, Real>::Zero();

    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            result(i, j) = a[i] * b[j];

    return result;
}




// Convenient aliases.
template <typename Real = float>
using Mat2 = MatN<2, 2, Real>;

template <typename Real = float>
using Mat3 = MatN<3, 3, Real>;

template <typename Real = float>
using Mat4 = MatN<4, 4, Real>;

template <typename Real = float>
using Mat6 = MatN<6, 6, Real>;
