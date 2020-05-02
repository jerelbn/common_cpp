#pragma once

namespace common
{

template <typename T, int R, int C, bool RowMajor = true>
class Matrix
{

public:
    Matrix() {}
    ~Matrix() {}

    Matrix(T *ptr)
    {
        memcpy(arr, ptr, sizeof(arr));
    }

    void operator=(Matrix<T, R, C, RowMajor> &m2)
    {
        for (int i = 0; i < R; ++i)
        {
            for (int j = 0; j < C; ++j)
            {
                (*this)(i, j) = m2(i, j);
            }
        }
    }

    T &operator()(int i, int j)
    {
        return arr[i * C + j];
    }

    template <int C2>
    Matrix<T, R, C2, RowMajor> operator*(Matrix<T, C, C2, RowMajor> &m2)
    {
        Matrix<T, R, C2, RowMajor> m;
        m.setZero();
        for (int i = 0; i < R; ++i)
        {
            for (int j = 0; j < C2; ++j)
            {
                for (int k = 0; k < C; ++k)
                {
                    m(i, j) += (*this)(i, k) * m2(k, j);
                }
            }
        }
        return m;
    }

    friend std::ostream &operator<<(std::ostream &os, Matrix<T, R, C, RowMajor> &m)
    {
        for (int i = 0; i < R; ++i)
        {
            for (int j = 0; j < C; ++j)
            {
                printf("%12.6f ", m(i, j));
            }
            printf("\n");
        }
        return os;
    }

    template <int R2, int C2>
    Matrix<T, R2, C2, RowMajor> block(int row, int col)
    {
        Matrix<T, R2, C2, RowMajor> m;
        for (int i = 0; i < R2; ++i)
        {
            for (int j = 0; j < C2; ++j)
            {
                m(i, j) = (*this)(row + i, col + j);
            }
        }
        return m;
    }

    void setZero()
    {
        memset(arr, T(0.0), sizeof(arr));
    }

    void setRandom(T scalar = T(1.0))
    {
        for (int i = 0; i < R; ++i)
        {
            for (int j = 0; j < C; ++j)
            {
                arr(i, j) = scalar * T(2.0) * T(rand()) / RAND_MAX - T(1.0);
            }
        }
    }

    static Matrix<T, R, C, RowMajor> zeros()
    {
        Matrix<T, R, C, RowMajor> m;
        memset(m.data(), T(0.0), m.size());
        return m;
    }

    static Matrix<T, R, C, RowMajor> random(T scalar = T(1.0))
    {
        Matrix<T, R, C, RowMajor> m;
        for (int i = 0; i < R; ++i)
        {
            for (int j = 0; j < C; ++j)
            {
                m(i, j) = scalar * T(2.0) * T(rand()) / RAND_MAX - T(1.0);
            }
        }
        return m;
    }

    const int rows() const { return R; }
    const int cols() const { return C; }
    const int size() const { return sizeof(arr); }
    T *const data() { return arr; }

private:
    T arr[R * C];
};

typedef Matrix<float, 1, 1> Matrix1f;
typedef Matrix<float, 2, 2> Matrix2f;
typedef Matrix<float, 3, 3> Matrix3f;
typedef Matrix<float, 4, 4> Matrix4f;
typedef Matrix<float, 5, 5> Matrix5f;
typedef Matrix<float, 6, 6> Matrix6f;
typedef Matrix<float, 7, 7> Matrix7f;
typedef Matrix<float, 8, 8> Matrix8f;
typedef Matrix<float, 9, 9> Matrix9f;

typedef Matrix<double, 1, 1> Matrix1d;
typedef Matrix<double, 2, 2> Matrix2d;
typedef Matrix<double, 3, 3> Matrix3d;
typedef Matrix<double, 4, 4> Matrix4d;
typedef Matrix<double, 5, 5> Matrix5d;
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 7, 7> Matrix7d;
typedef Matrix<double, 8, 8> Matrix8d;
typedef Matrix<double, 9, 9> Matrix9d;

typedef Matrix<float, 1, 1> Vector1f;
typedef Matrix<float, 2, 1> Vector2f;
typedef Matrix<float, 3, 1> Vector3f;
typedef Matrix<float, 4, 1> Vector4f;
typedef Matrix<float, 5, 1> Vector5f;
typedef Matrix<float, 6, 1> Vector6f;
typedef Matrix<float, 7, 1> Vector7f;
typedef Matrix<float, 8, 1> Vector8f;
typedef Matrix<float, 9, 1> Vector9f;

typedef Matrix<double, 1, 1> Vector1d;
typedef Matrix<double, 2, 1> Vector2d;
typedef Matrix<double, 3, 1> Vector3d;
typedef Matrix<double, 4, 1> Vector4d;
typedef Matrix<double, 5, 1> Vector5d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 7, 1> Vector7d;
typedef Matrix<double, 8, 1> Vector8d;
typedef Matrix<double, 9, 1> Vector9d;

} // namespace common