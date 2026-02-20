#include "Point.h"
#include "Interval.h"
#include "IOptProblem.hpp"

template<class T, std::size_t N>
std::vector<T> normalize(const TPoint<T, N>& point, const TPoint<T, N>& lowerBound, const TPoint<T, N>& upperBound)
{
  TPoint<T, N> normalized;
  for (std::size_t i = 0; i < N; ++i) {
    T range = upperBound[i] - lowerBound[i];
    if (std::abs(range) < 1e-10) {
      throw std::domain_error("Zero range in dimension");
    }
    normalized[i] = (point[i] - lowerBound[i]) / range;
  }
  return normalized.toVector();
}

template<class T, std::size_t N>
std::vector<T> denormalize(const TPoint<T, N>& point, const TPoint<T, N>& lowerBound, const TPoint<T, N>& upperBound)
{
  TPoint<T, N> denormalized;
  for (std::size_t i = 0; i < N; ++i) {
    T range = upperBound[i] - lowerBound[i];
    denormalized[i] = lowerBound[i] + point[i] * range;
  }
  return denormalized.toVector();
}

template<class T, size_t N>
class TAlgorithm
{
private:
    IOptProblem* func;
    TInterval<double> interval;
    double eps;
    size_t r;
    double R;
    double M;
    double L;
    double zi;
    double zi1;
    double xi;
    double xi1;
    TPoint<T, N> lowerBound;
    TPoint<T, N> upperBound;
    size_t indexInteravlWhithMaxR;
    size_t indexInteravlWhithMinR;
    double newPoint;
    size_t iteration;
public:
    TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
              size_t r_, IOptProblem* func_);
    ~TAlgorithm();
    
    // Единый метод для решения. По умолчанию ищет минимум.
    std::vector<T> Solve(bool isMinimize = true); 
};

template <class T, size_t N>
inline TAlgorithm<T, N>::TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
              size_t r_, IOptProblem* func_)
{
    if(eps >= 1)
        throw std::invalid_argument("eps >= 1");
    if(r_ < 2)
        throw std::invalid_argument("r < 2");
    eps = eps_;
    r = r_;
    func = func_;
    iteration = 0;
    lowerBound = lowerBound_;
    upperBound = upperBound_;
}

template <class T, size_t N>
inline TAlgorithm<T, N>::~TAlgorithm()
{
}

template <class T, size_t N>
inline std::vector<T> TAlgorithm<T, N>::Solve(bool isMinimize)
{
    // Определяем множитель. Для минимума 1.0, для максимума -1.0
    T sign = isMinimize ? 1.0 : -1.0;

    T resX = denormalize(TPoint<T,N>(interval.getLeft(0)), lowerBound, upperBound)[0];
    
    // Умножаем результат функции на sign
    T resZ_internal = sign * func->ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(0)), lowerBound, upperBound));

    do
    {
        M = 0;
        iteration++;
        
        // 1. Вычисление константы Липшица 
        for (size_t i = 0; i < interval.size(); i++)
        {
            zi = sign * func->ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(i)), lowerBound, upperBound));
            zi1 = sign * func->ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(i)), lowerBound, upperBound));
            xi = interval.getRight(i);
            xi1 = interval.getLeft(i);
            
            M = std::max(M, std::abs(zi - zi1) / (xi - xi1));
            
            if (zi < resZ_internal)
            {
                resX = denormalize(TPoint<T,N>(xi), lowerBound, upperBound)[0];
                resZ_internal = zi;
            }
        }
        
        if (M == 0)
            L = 1;
        else
            L = r * M;
            
        // 2. Вычисление характеристик R для каждого интервала
        for (size_t i = 0; i < interval.size(); i++)
        {
            zi = sign * func->ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(i)), lowerBound, upperBound));
            zi1 = sign * func->ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(i)), lowerBound, upperBound));
            xi = interval.getRight(i);
            xi1 = interval.getLeft(i);
            
            R = L * (xi - xi1) + ((zi - zi1) * (zi - zi1) / (L * (xi - xi1))) - 2 * (zi + zi1);
            interval.setIntervalR(i, R);
        }
        
        indexInteravlWhithMaxR = interval.getMaxRIntervalIndex();
        
        // 3. Выбор точки для нового испытания
        zi = sign * func->ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(indexInteravlWhithMaxR)), lowerBound, upperBound));
        zi1 = sign * func->ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(indexInteravlWhithMaxR)), lowerBound, upperBound));
        xi = interval.getRight(indexInteravlWhithMaxR);
        xi1 = interval.getLeft(indexInteravlWhithMaxR);
        
        newPoint = ((xi1 + xi) / 2) - ((zi - zi1) / (2 * L));
        T newZ_internal = sign * func->ComputeFunction(denormalize(TPoint<T,N>(newPoint), lowerBound, upperBound));
        
        if (newZ_internal < resZ_internal)
        {
            resX = denormalize(TPoint<T,N>(newPoint), lowerBound, upperBound)[0];
            resZ_internal = newZ_internal;
        }
        
        interval.split(newPoint);
    }
    while ((interval.getLength(interval.getMaxRIntervalIndex()) > eps) && (iteration < 1000));

    // 4. Финальная проверка на последнем разбитом интервале
    size_t maxIndex = interval.getMaxRIntervalIndex();
    double x0 = interval.getLeft(maxIndex);
    double x1 = interval.getRight(maxIndex);
    T tmp1 = denormalize(TPoint<T, N>(x0), lowerBound, upperBound)[0];
    T tmp2 = denormalize(TPoint<T, N>(x1), lowerBound, upperBound)[0];
    T tmpX = (tmp1 + tmp2) / 2;
    std::vector<T> xd(1);
    xd[0] = tmpX;
    
    T tmpZ_internal = sign * func->ComputeFunction(xd);
    if (tmpZ_internal < resZ_internal)
    {
        resX = tmpX;
        resZ_internal = tmpZ_internal;
    }
    
    // Возвращаем знак функции в исходное состояние, если искали максимум
    T finalZ = isMinimize ? resZ_internal : -resZ_internal;
    
    return std::vector<T>({finalZ, resX, (double)iteration});
}