#include <array>
#include <iostream>
#pragma once

/**
* xType - тип аргумента x.
* yType - тип значения функции y
* N - количество точек для интерполяции
*
* Рекомедую обратить внимание. Разность (xType - xType) не обязана быть типом xType
*/
template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator
{

    std::array<xType, N> func_points;
    std::array<yType, N> diff;

    public:
    NewtonInterpolator(const std::array<xType, N> &points, const std::array<yType, N>& values) noexcept
    : func_points{points}, diff(values)
    {

        for (std::size_t i = 0; i < N - 1; i++)
        {
            for (std::size_t j = N - 1; j > i; j--)
            {
                diff[j] = (diff[j] - diff[j - 1]) / (points[j] - points[j - 1 - i]);
            }
        }
    };

    yType interpolate(const xType& x) const noexcept
    {
        yType res = diff[N-1];
        for (std::size_t i = N - 1; i > 0; i--)
        {
            res = std::fma(res, (x - func_points[i - 1]), diff[i - 1]);
        }
        return res;
    }
};