#pragma once
#include <array>
#include<eigen3/Eigen/Dense>
#include <iostream>

template<typename RealType, unsigned int N>
struct DerivativeCoef 
{
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};

template<typename RealType, unsigned int N>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept
{
    Eigen::Matrix<RealType, N, N> A;
    Eigen::Vector<RealType, N> b;
    Eigen::Vector<RealType, N> x;

    for (std::size_t i = 0; i < N; i++)                         //заполняем матрицу коэфф без первой строки и столбца, центральный коэфф считаем отдельно далее
    {
        A(0, i) = points[i];
        b(i) = 0;                                                           //для поиска первой производной
    }

    b(0) = 1; 
    for (std::size_t i = 1; i < N; i++)
    {
        for (std::size_t j = 0; j < N; j++)
        {
            A(i, j) = A(i - 1, j) * points[j];
        }
    }

    x = A.householderQr().solve(b);
    RealType central = -x.sum();

    std::array<RealType, N> otherx;
    for (std::size_t i = 0; i < N; i++)
    {
        otherx[i] = x(i);
    }

    return {central, otherx};
};
