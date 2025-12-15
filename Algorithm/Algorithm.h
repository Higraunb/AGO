#include "Point.h"
#include "Interval.h"
#include "IOptProblem.hpp"

template<class T, std::size_t N>
std::vector<T>& normalize(const TPoint<T, N>& point)
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
std::vector<T>& denormalize(const TPoint<T, N>& point)
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
    size_t indexInteravlWhithMaxR;
    size_t indexInteravlWhithMinR;
    double newPoint;
    size_t iteration;
public:
    TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
              size_t r_, IOptProblem* func_);
    ~TAlgorithm();
    void AGPStronginaMin();
    void AGPStronginaMax();
};
template <class T, size_t N>
inline TAlgorithm<T, N>::TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
              size_t r_, IOptProblem* func_)
{
    // if(lowerBound_ >= upperBound_)
    //     throw std::invalid_argument("lowerBound >= upperBound");
    if(eps_ >= 1)
        throw std::invalid_argument("eps >= 1");
    if(r_ < 2)
        throw std::invalid_argument("r < 2");
    eps = eps_;
    r = r_;
    func = func_;
    iteration = 0;
}

template <class T, size_t N>
inline TAlgorithm<T, N>::~TAlgorithm()
{
}

template <class T, size_t N>
inline void TAlgorithm<T, N>::AGPStronginaMin()
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
        zi = func.ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(i))));
        zi1 = func.ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(i))));
        xi = interval.getRigth(i);
        xi1 = interval.getLeft(i);
        M = max(M, abs(zi - zi1) / (xi - xi1));
    }
    if (M == 0)
        L = 1;
    else
        L = r * M;
    for (size_t i = 0; i < interval.size(); i++)
    {
        //zi = func.evaluateNormalized(TPoint<T,N>(interval.getRigth(i)));
        //zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(i)));
        zi = func.ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(i))));
        zi1 = func.ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(i))));
        xi = interval.getRigth(i);
        xi1 = interval.getLeft(i);
        R = L * (xi - xi1) + (pow(zi - zi1, 2) / (L * (xi - xi1))) - 2 * (zi + zi1);
        interval.setIntervalR(i, R);
    }
    indexInteravlWhithMaxR = interval.getMaxRIntervalIndex();
    //zi = func.evaluateNormalized(TPoint<T,N>(interval.getRigth(indexInteravlWhithMaxR)));
    //zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(indexInteravlWhithMaxR)));
    zi = func.ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(indexInteravlWhithMaxR))));
    zi1 = func.ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(indexInteravlWhithMaxR))));
    xi = interval.getRigth(indexInteravlWhithMaxR);
    xi1 = interval.getLeft(indexInteravlWhithMaxR);
    newPoint = ((xi1 + xi) / 2) - ((zi - zi1) / (2 * L));
    interval.split(newPoint);
    } 
    while (interval.getLength(interval.getMaxRIntervalIndex()) > eps);
    size_t maxIndex = interval.getMaxRIntervalIndex();
    double x0 = interval.getLeft(maxIndex);
    double x1 = interval.getRigth(maxIndex);
    TPoint<T, N> resX = (denormalize(TPoint<T, N>(x0))[0] + denormalize(TPoint<T, N>(x1))[0])/2;
    T resY = func.ComputeFunction(resX.data());
    cout << "xr = " << resX << endl;
    cout << "fr = " << resY << endl;
    cout << "Iterations = " << iteration << endl;
}

template <class T, size_t N>
inline void TAlgorithm<T, N>::AGPStronginaMax()
{
    //func.setMaximization(true);
    do
    {
    M = 0;
    iteration++;
    for (size_t i = 0; i < interval.size(); i++)
    {
        //zi = func.evaluateNormalized(TPoint<T,N>(interval.getRigth(i)));
        //zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(i)));
        zi = -func.ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(i))));
        zi1 = -func.ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(i))));
        xi = interval.getRigth(i);
        xi1 = interval.getLeft(i);
        M = max(M, abs(zi - zi1) / (xi - xi1));
    }
    if (M == 0)
        L = 1;
    else
        L = r * M;
    for (size_t i = 0; i < interval.size(); i++)
    {
        //zi = func.evaluateNormalized(TPoint<T,N>(interval.getRigth(i)));
        //zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(i)));
        zi = -func.ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(i))));
        zi1 = -func.ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(i))));
        xi = interval.getRigth(i);
        xi1 = interval.getLeft(i);
        R = L * (xi - xi1) + (pow(zi - zi1, 2) / (L * (xi - xi1))) - 2 * (zi + zi1);
        interval.setIntervalR(i, R);
    }
    indexInteravlWhithMaxR = interval.getMaxRIntervalIndex();
    //zi = func.evaluateNormalized(TPoint<T,N>(interval.getRigth(indexInteravlWhithMaxR)));
    //zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(indexInteravlWhithMaxR)));
    zi = -func.ComputeFunction(denormalize(TPoint<T,N>(interval.getRight(indexInteravlWhithMaxR))));
    zi1 = -func.ComputeFunction(denormalize(TPoint<T,N>(interval.getLeft(indexInteravlWhithMaxR))));
    xi = interval.getRigth(indexInteravlWhithMaxR);
    xi1 = interval.getLeft(indexInteravlWhithMaxR);
    newPoint = ((xi1 + xi) / 2) - ((zi - zi1) / (2 * L));
    interval.split(newPoint);
    } 
    while (interval.getLength(interval.getMaxRIntervalIndex()) > eps);
    size_t maxIndex = interval.getMaxRIntervalIndex();
    double x0 = interval.getLeft(maxIndex);
    double x1 = interval.getRigth(maxIndex);
    TPoint<T, N> resX = (denormalize(TPoint<T, N>(x0))[0] + denormalize(TPoint<T, N>(x1))[0])/2;
    T resY = func.ComputeFunction(resX.data());
    cout << "xr = " << resX << endl;
    cout << "fr = " << -1 * resY << endl;
    cout << "Iterations = " << iteration << endl;
}
