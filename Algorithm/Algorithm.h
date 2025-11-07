
#include "Function.h"
#include <math.h>
template<class T, size_t N>
class TAlgorithm
{
private:
    using FunctionType = std::function<T(const TPoint<T, N>&)>;
    TFunction<T, N> func;
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
              size_t r_, FunctionType func_);
    TAlgorithm(double eps_, size_t r_, TFunction<T, N> func_);
    ~TAlgorithm();
    double AGPStronginaMin();
    void AGPStronginaMax();
};

template <class T, size_t N>
inline TAlgorithm<T, N>::TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
              size_t r_, FunctionType func_)
{
    // if(lowerBound_ >= upperBound_)
    //     throw std::invalid_argument("lowerBound >= upperBound");
    if(eps_ >= 1)
        throw std::invalid_argument("eps >= 1");
    if(r_ < 2)
        throw std::invalid_argument("r < 2");
    eps = eps_;
    r = r_;
    func = TFunction<T,N>(func_, lowerBound_, upperBound_);
    iteration = 0;
}

template <class T, size_t N>
inline TAlgorithm<T, N>::TAlgorithm(double eps_, size_t r_, TFunction<T, N> func_)
{
    if(eps_ >= 1)
        throw std::invalid_argument("eps >= 1");
    if(r_ < 2)
        throw std::invalid_argument("r < 2");
    eps = eps_;
    r = r_;
    func = func_;
}

template <class T, size_t N>
inline TAlgorithm<T, N>::~TAlgorithm()
{
}

template <class T, size_t N>
inline double TAlgorithm<T, N>::AGPStronginaMin()
{
    func.setMaximization(false);
    do
    {
    M = 0;
    iteration++;
    for (size_t i = 0; i < interval.size(); i++)
    {
        zi = func.evaluateNormalized(TPoint<T,N>(interval.getRigth(i)));
        zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(i)));
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
        zi = func.evaluateNormalized(TPoint<T,N>(interval.getRigth(i)));
        zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(i)));
        xi = interval.getRigth(i);
        xi1 = interval.getLeft(i);
        R = L * (xi - xi1) + (pow(zi - zi1, 2) / (L * (xi - xi1))) - 2 * (zi + zi1);
        interval.setIntervalR(i, R);
    }
    indexInteravlWhithMaxR = interval.getMaxRIntervalIndex();
    zi = func.evaluateNormalized(TPoint<T,N>(interval.getRigth(indexInteravlWhithMaxR)));
    zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(indexInteravlWhithMaxR)));
    xi = interval.getRigth(indexInteravlWhithMaxR);
    xi1 = interval.getLeft(indexInteravlWhithMaxR);
    newPoint = ((xi1 + xi) / 2) - ((zi - zi1) / (2 * L));
    interval.split(newPoint);
    } 
    while (interval.getLength(interval.getMaxRIntervalIndex()) > eps);
    size_t maxIndex = interval.getMaxRIntervalIndex();
    double x0 = interval.getLeft(maxIndex);
    double x1 = interval.getRigth(maxIndex);
    TPoint<T, N> resX = (func.denormalize(TPoint<T, N>(x0)) + func.denormalize(TPoint<T, N>(x1)))/2;
    double resY = func(resX);
    cout << "xr = " << resX << endl;
    cout << "fr = " << resY << endl;
    cout << "Iterations = " << iteration << endl;
    return resY;
}

template <class T, size_t N>
inline void TAlgorithm<T, N>::AGPStronginaMax()
{
    func.setMaximization(true);
    do
    {
    M = 0;
    iteration++;
    for (size_t i = 0; i < interval.size(); i++)
    {
        zi = func.evaluateNormalized(TPoint<T,N>(interval.getRigth(i)));
        zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(i)));
        xi = interval.getRigth(i);
        xi1 = interval.getLeft(i);
        M = std::max(M, abs(zi - zi1) / (xi - xi1));
    }
    if (M == 0)
        L = 1;
    else
        L = r * M;
    for (size_t i = 0; i < interval.size(); i++)
    {
        zi = func.evaluateNormalized(TPoint<T,N>(interval.getRigth(i)));
        zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(i)));
        xi = interval.getRigth(i);
        xi1 = interval.getLeft(i);
        R = L * (xi - xi1) + (pow(zi - zi1, 2) / (L * (xi - xi1))) - 2 * (zi + zi1);
        interval.setIntervalR(i, R);
    }
    indexInteravlWhithMaxR = interval.getMaxRIntervalIndex();
    zi = func.evaluateNormalized(TPoint<T,N>(interval.getRigth(indexInteravlWhithMaxR)));
    zi1 = func.evaluateNormalized(TPoint<T,N>(interval.getLeft(indexInteravlWhithMaxR)));
    xi = interval.getRigth(indexInteravlWhithMaxR);
    xi1 = interval.getLeft(indexInteravlWhithMaxR);
    newPoint = ((xi1 + xi) / 2) - ((zi - zi1) / (2 * L));
    interval.split(newPoint);
    } 
    while (interval.getLength(interval.getMaxRIntervalIndex()) > eps);
    size_t maxIndex = interval.getMaxRIntervalIndex();
    double x0 = interval.getLeft(maxIndex);
    double x1 = interval.getRigth(maxIndex);
    TPoint<T, N> resX = (func.denormalize(TPoint<T, N>(x0)) + func.denormalize(TPoint<T, N>(x1)))/2;
    T resY = func(resX);
    std::cout << "xr = " << resX << endl;
    std::cout << "fr = " << -1 * resY << endl;
    std::cout << "Iterations = " << iteration << endl;
}