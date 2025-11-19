#include "Function.h"
#include <map>
using namespace std;

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
    std::vector<double> allPoint;
    TPoint<T, N> resX;
    T resY;
public:
    TAlgorithm();
    TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
              size_t r_, FunctionType func_);
    ~TAlgorithm();

    T GetLowerBound();
    T GetUpperBound();
    FunctionType GetFunc();
    string GetInfo();

    vector<double> GetAllPoint();
    void AGPStronginaMin();
    void AGPStronginaMax();
};

template<class T, size_t N>
TAlgorithm<T, N>::TAlgorithm()
{
    func = TFunction<T,N>([](const TPoint<T, N>& x){return x[0];},TPoint<T,N>(0),TPoint<T,N>(1));
    eps = 0.01;
    r = 2;
    iteration = 0;
}

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
inline TAlgorithm<T, N>::~TAlgorithm()
{
}

template<class T, size_t N>
inline string TAlgorithm<T, N>::GetInfo()
{
    string info = "xr = " + to_string(resX[0]) + '\n' + "fr = "
                  + to_string(resY) + '\n' + "Iterations = "
                  + to_string(iteration);
    return info;
}

template<class T, size_t N>
inline vector<double> TAlgorithm<T, N>::GetAllPoint()
{
    return allPoint;
}

template <class T, size_t N>
inline void TAlgorithm<T, N>::AGPStronginaMin()
{
    func.setMaximization(false);
    allPoint.clear();
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
    TPoint<T,N> realNewPoint = func.denormalize(TPoint<T,N>(newPoint));
    allPoint.push_back(realNewPoint[0]);
    interval.split(newPoint);
    } 
    while (interval.getLength(interval.getMaxRIntervalIndex()) > eps);

    size_t maxIndex = interval.getMaxRIntervalIndex();
    double x0 = interval.getLeft(maxIndex);
    double x1 = interval.getRigth(maxIndex);
    resX = (func.denormalize(TPoint<T, N>(x0)) + func.denormalize(TPoint<T, N>(x1)))/2;
    resY = func(resX);
    // cout << "xr = " << resX << endl;
    // cout << "fr = " << resY << endl;
    // cout << "Iterations = " << iteration << endl;
}

template <class T, size_t N>
inline void TAlgorithm<T, N>::AGPStronginaMax()
{
    allPoint.clear();
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
    TPoint<T,N> realNewPoint = func.denormalize(TPoint<T,N>(newPoint));
    allPoint.push_back(realNewPoint[0]);
    interval.split(newPoint);
    } 
    while (interval.getLength(interval.getMaxRIntervalIndex()) > eps);

    size_t maxIndex = interval.getMaxRIntervalIndex();
    double x0 = interval.getLeft(maxIndex);
    double x1 = interval.getRigth(maxIndex);
    resX = (func.denormalize(TPoint<T, N>(x0)) + func.denormalize(TPoint<T, N>(x1)))/2;
    resY = func(resX);
    // cout << "xr = " << resX << endl;
    // cout << "fr = " << -1 * resY << endl;
    // cout << "Iterations = " << iteration << endl;
}
