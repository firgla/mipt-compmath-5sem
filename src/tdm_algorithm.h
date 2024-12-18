#include <iostream>
#include <vector>

struct Tridiagonal_Matrix
{
    std::vector<double> A;
    std::vector<double> B;
    std::vector<double> C;

};

std::vector<double> Algorithm(Tridiagonal_Matrix& matrix, std::vector<double>& f);
