#include "tdm_algorithm.h"

std::vector<double> Algorithm(Tridiagonal_Matrix& matrix, std::vector<double>& f)
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
}
