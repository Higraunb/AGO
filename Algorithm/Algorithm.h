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
  return normalized.data();
}

template<class T, std::size_t N>
std::vector<T> denormalize(const TPoint<T, N>& point, const TPoint<T, N>& lowerBound, const TPoint<T, N>& upperBound)
{
  TPoint<T, N> denormalized;
  for (std::size_t i = 0; i < N; ++i) {
    T range = upperBound[i] - lowerBound[i];
    denormalized[i] = lowerBound[i] + point[i] * range;
  }
  return denormalized.data();
}

template<class T, size_t N>
class TAlgorithm
{
private:
    //using FunctionType = std::function<T(const TPoint<T, N>&)>;
    //TFunction<T, N> func;
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
    std::vector<T> AGPStronginaMin();
    std::vector<T> AGPStronginaMax();
};
template <class T, size_t N>
inline TAlgorithm<T, N>::TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
              size_t r_, IOptProblem* func_)
{
    if(eps_ >= 1)
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
inline std::vector<T> TAlgorithm<T, N>::AGPStronginaMin()
{
    //func.setMaximization(false);
    do
    {
    M = 0;
    iteration++;
    for (size_t i = 0; i < interval.size(); i++)
    {
        //zi = func.evaluateNormalized(TPoint<T,N>(interval.getRight(i)));
        //zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(i)));
        zi = func->ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(i)), lowerBound, upperBound));
        zi1 = func->ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(i)), lowerBound, upperBound));
        xi = interval.getRight(i);
        xi1 = interval.getLeft(i);
        M = std::max(M, abs(zi - zi1) / (xi - xi1));
    }
    if (M == 0)
        L = 1;
    else
        L = r * M;
    for (size_t i = 0; i < interval.size(); i++)
    {
        //zi = func.evaluateNormalized(TPoint<T,N>(interval.getRight(i)));
        //zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(i)));
        zi = func->ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(i)), lowerBound, upperBound));
        zi1 = func->ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(i)), lowerBound, upperBound));
        xi = interval.getRight(i);
        xi1 = interval.getLeft(i);
        R = L * (xi - xi1) + (pow(zi - zi1, 2) / (L * (xi - xi1))) - 2 * (zi + zi1);
        interval.setIntervalR(i, R);
    }
    indexInteravlWhithMaxR = interval.getMaxRIntervalIndex();
    //zi = func.evaluateNormalized(TPoint<T,N>(interval.getRight(indexInteravlWhithMaxR)));
    //zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(indexInteravlWhithMaxR)));
    zi = func->ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(indexInteravlWhithMaxR)), lowerBound, upperBound));
    zi1 = func->ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(indexInteravlWhithMaxR)), lowerBound, upperBound));
    xi = interval.getRight(indexInteravlWhithMaxR);
    xi1 = interval.getLeft(indexInteravlWhithMaxR);
    newPoint = ((xi1 + xi) / 2) - ((zi - zi1) / (2 * L));
    interval.split(newPoint);
    } 
    while (interval.getLength(interval.getMaxRIntervalIndex()) > eps);
    size_t maxIndex = interval.getMaxRIntervalIndex();
    double x0 = interval.getLeft(maxIndex);
    double x1 = interval.getRight(maxIndex);
    T tmp1 = denormalize(TPoint<T, N>(x0), lowerBound, upperBound)[0];
    T tmp2 = denormalize(TPoint<T, N>(x1), lowerBound, upperBound)[0];
    T resX = (tmp1 + tmp2)/2;
    std::vector<T> xd(1);
    xd[0] = resX;
    T resY = func->ComputeFunction(xd);
    return std::vector<T>({resY, resX, (double) iteration});
    // std::cout << "xr = " << resX << std::endl;
    // std::cout << "fr = " << resY << std::endl;
    // std::cout << "Iterations = " << iteration << std::endl;
}

template <class T, size_t N>
inline std::vector<T> TAlgorithm<T, N>::AGPStronginaMax()
{
    //func.setMaximization(true);
    do
    {
    M = 0;
    iteration++;
    for (size_t i = 0; i < interval.size(); i++)
    {
        //zi = func.evaluateNormalized(TPoint<T,N>(interval.getRight(i)));
        //zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(i)));
        zi = -func->ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(i)), lowerBound, upperBound));
        zi1 = -func->ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(i)), lowerBound, upperBound));
        xi = interval.getRight(i);
        xi1 = interval.getLeft(i);
        M = std::max(M, abs(zi - zi1) / (xi - xi1));
    }
    if (M == 0)
        L = 1;
    else
        L = r * M;
    for (size_t i = 0; i < interval.size(); i++)
    {
        //zi = func.evaluateNormalized(TPoint<T,N>(interval.getRight(i)));
        //zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(i)));
        zi = -func->ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(i)), lowerBound, upperBound));
        zi1 = -func->ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(i)), lowerBound, upperBound));
        xi = interval.getRight(i);
        xi1 = interval.getLeft(i);
        R = L * (xi - xi1) + (pow(zi - zi1, 2) / (L * (xi - xi1))) - 2 * (zi + zi1);
        interval.setIntervalR(i, R);
    }
    indexInteravlWhithMaxR = interval.getMaxRIntervalIndex();
    //zi = func.evaluateNormalized(TPoint<T,N>(interval.getRight(indexInteravlWhithMaxR)));
    //zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(indexInteravlWhithMaxR)));
    zi = -func->ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(indexInteravlWhithMaxR)), lowerBound, upperBound));
    zi1 = -func->ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(indexInteravlWhithMaxR)), lowerBound, upperBound));
    xi = interval.getRight(indexInteravlWhithMaxR);
    xi1 = interval.getLeft(indexInteravlWhithMaxR);
    newPoint = ((xi1 + xi) / 2) - ((zi - zi1) / (2 * L));
    interval.split(newPoint);
    } 
    while (interval.getLength(interval.getMaxRIntervalIndex()) > eps);
    size_t maxIndex = interval.getMaxRIntervalIndex();
    double x0 = interval.getLeft(maxIndex);
    double x1 = interval.getRight(maxIndex);
    T resX = (denormalize(TPoint<T, N>(x0), lowerBound, upperBound)[0] + denormalize(TPoint<T, N>(x1), lowerBound, upperBound)[0])/2;
    std::vector<T> xd(1);
    xd[0] = resX;
    T resY = func->ComputeFunction(xd);
    return std::vector({resY, resX, (double) iteration});
    // std::cout << "xr = " << resX << std::endl;
    // std::cout << "fr = " << -1 * resY << std::endl;
    // std::cout << "Iterations = " << iteration << std::endl;
}
