#include <vector>
#include <type_traits>

/** класс для работы с трехдиагональной матрицей **/
template<typename Type>
class ThreeDiagonalMatrix
{
    public:
    std::vector<double> A;
    std::vector<double> B;
    std::vector<double> C;

    ThreeDiagonalMatrix(
        const std::vector<Type> &a,
        const std::vector<Type> &b,
        const std::vector<Type> &c) : A(a), B(b), C(c)
        {}

};


template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

/** Функция для решения методм  прогонки **/
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(ThreeDiagonalMatrix<mType>& matrix, std::vector<cType>& f)
{
    std::vector<double> q = {0.0};
    std::vector<double> p = {0.0};
    std::vector<double> x;

    matrix.A.insert(matrix.A.begin(), 0.0);
    matrix.C.insert(matrix.C.end(), 0.0);

    for(int i = 0; i < matrix.B.size() - 1; i++)
    {
        double q_i1 = (f[i] - matrix.A[i] * q[i]) / (matrix.A[i] * p[i] + matrix.B[i]);
        q.emplace_back(q_i1);

        double p_i1 = (- matrix.C[i]) / (matrix.A[i] * p[i] + matrix.B[i]);
        p.emplace_back(p_i1);

    }

    double x_n = (f.back() - matrix.A.back() * q.back()) / (matrix.A.back() * p.back() + matrix.B.back());
    x.emplace_back(x_n);

    for (int i = matrix.B.size() - 1; i > 0; i--)
    {
        double x_i1 = p[i] * x.front() + q[i];
        x.insert(x.begin(), x_i1);
    }

    return x;
};

template<typename Type>
struct SplinePolynomCoeffs 
{
    Type a, b, c, d;
};

/**
* xType - тип аргумента x.
* yType - тип значения функции y
*/

template<typename xType, typename yType>
class CubicSpline 
{
    std::vector<xType> func_points;
    std::vector<SplinePolynomCoeffs<DivisType<xType, yType>>> func_coeffs;
    
    public:
    CubicSpline(const std::vector<xType> &points, const std::vector<yType>& values)
    :func_points(points)
    {
        std::size_t n = func_points.size();
        func_coeffs.resize(n - 1);
        std::vector<xType> h(n - 1);
        std::vector<yType> u(n - 1);
        std::vector<DivisType<xType, yType>> column(n - 2);
        
        for(std::size_t i = 0; i < n - 1; i++)
        {
            h[i] = func_points[i + 1] - func_points[i];
            u[i] = (values[i + 1] - values[i]) / h[i];
        }
        
        for (std::size_t i = 0; i < n - 2; i++) 
        {
            column[i] = 6 * (u[i + 1] - u[i]) / (h[i + 1] + h[i]);
        }

        std::vector<xType> b(n - 2, 2);
        std::vector<xType> a(n - 3), c(n - 3);

        for (std::size_t i = 0; i < n - 3; i++) 
        {
            c[i] = h[i + 1] / (h[i + 1] + h[i]);
            a[i] = h[i] / (h[i + 1] + h[i]);
        }
       
        ThreeDiagonalMatrix<xType> matrix(a, b, c);
        std::vector<DivisType<xType, yType>> C_coeffs = solve(matrix, column);

        for (std::size_t i = 0; i < n - 1; i++) 
        {
            func_coeffs[i].a = values[i + 1];
        }

        func_coeffs[n - 2].c = 0;
        for (std::size_t i = 0; i < n - 2; i++)
        {
            func_coeffs[i].c = C_coeffs[i];
        }

      
        func_coeffs[0].d = func_coeffs[0].c / h[0];

        for (std::size_t i = 1; i < n-1; i++) 
        {
            func_coeffs[i].d = (func_coeffs[i].c - func_coeffs[i - 1].c) / h[i];
        }

        for (std::size_t i = 0; i < n - 1; i++) 
        {
            func_coeffs[i].b = func_coeffs[i].c / 2 * h[i] - func_coeffs[i].d/6 * h[i] * h[i]+ u[i];
        }
    }
                        
    yType interpolate(const xType& x) const noexcept
    {
        for (std::size_t i = 0; i < func_coeffs.size(); ++i) 
        {
            if (func_points[i] <= x && x <= func_points[i + 1]) 
                return func_coeffs[i].a + (x - func_points[i + 1]) * (func_coeffs[i].b  +  (x - func_points[i + 1]) * (func_coeffs[i].c / 2 + (x - func_points[i + 1]) * func_coeffs[i].d / 6));
        }        
    }

};

