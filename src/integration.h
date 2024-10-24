#include <array>
#include <type_traits>

template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> 
{
    using Argument = Arg;
};

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());

template <std::size_t N>
struct GaussianKvadrature
{
    static constexpr std::array<double, N> points, weights;
};

template <>
struct GaussianKvadrature<2>
{
    static constexpr std::array<double, 2> points = {-0.57735, 0.57735} , weights = {1, 1};
};

template <>
struct GaussianKvadrature<3>
{
    static constexpr std::array<double, 3> points = {-0.77459666924148338, 0, 0.77459666924148338}, weights = {5.0 / 9, 8.0 / 9, 5.0 / 9};
};

template <>
struct GaussianKvadrature<4>
{
    static constexpr std::array<double, 4> points = {-0.86113631159405257254, -0.33998104358485625731, 0.33998104358485625731, 0.86113631159405257254};
    static constexpr std::array<double, 4> weights = {0.34785484513745390522, 0.65214515486254631682, 0.65214515486254631682, 0.34785484513745390522};
};

/* Функция производит интегрирование на одном отрезке */
template<typename Callable, std::size_t N>
decltype(auto) integrate(const Callable& func, const typename ArgumentGetter<Callable>::Argument& start, const typename ArgumentGetter<Callable>::Argument& end)
{
    typename ArgumentGetter<Callable>::Argument answer = GaussianKvadrature<N>::weights[0] * func((end + start) / 2 + GaussianKvadrature<N>::points[0] * (end - start) / 2);

    for (size_t i = 1; i < N; i++)
    {
        answer += GaussianKvadrature<N>::weights[i] * func((end + start) / 2 + GaussianKvadrature<N>::points[i] * (end - start) / 2);
    }

    return answer * (end - start) / 2;
};

/* Функция производит интегрирование, разбивая отрезок на подотрезки длиной не более dx */
template<typename Callable, std::size_t N>
decltype(auto) integrate(const Callable& func, const typename ArgumentGetter<Callable>::Argument& start, const typename ArgumentGetter<Callable>::Argument& end, const Dif<typename ArgumentGetter<Callable>::Argument>& dx)
{
    typename ArgumentGetter<Callable>::Argument answer = integrate<Callable, N>(func, start, start + dx);
    typename ArgumentGetter<Callable>::Argument curr_position = start + dx;
    while (curr_position + dx < end)
    {
        answer += integrate<Callable,  N>(func, curr_position, curr_position + dx);
        curr_position += dx;
    }

    answer += integrate<Callable, N>(func, curr_position, end);

    return answer;
};
