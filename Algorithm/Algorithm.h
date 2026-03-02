#include "Point.h"
#include "Interval.h"
#include "IOptProblem.hpp"
#include "Evolvent.hpp"
#include <queue>

/*template<class T, std::size_t N>
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
}*/

template<class T, size_t N>
class TAlgorithm
{
private:
    IOptProblem* func;
    TInterval<double> interval;
    double eps;
    double r;
    double R;
    double M;
    double L;
    double zi;
    double zi1;
    double xi;
    double xi1;
    double newX;
    ags::Evolvent evolvent;
    TPoint<T, N> lowerBound;
    TPoint<T, N> upperBound;
    size_t indexInteravlWhithMaxR;
    size_t indexInteravlWhithMinR;
    size_t iteration;
    size_t CalculateR(const size_t& LeipshitzConstant, const double& xRight, const double& xLeft, const double& zRight, const double& zLeft);
    
public:
    TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
              double r_, IOptProblem* func_, int tightness_);
    
    ~TAlgorithm();
    
    std::vector<T> Solve(size_t maxInteration, bool isMinimize = true); 
};

template <class T, size_t N>
inline size_t TAlgorithm<T, N>::CalculateR(const size_t& LeipshitzConstant, const double& xRight, const double& xLeft, const double& zRight, const double& zLeft)
{
    double res = pow((xRight - xLeft), 1.0 / N) + ((zRight - zLeft) * (zRight - zLeft) / (L * pow((xRight - xLeft), 1.0 / N))) - 2 * (zRight + zLeft);
    return res;
}

template <class T, size_t N>
inline TAlgorithm<T, N>::TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
                                    double r_, IOptProblem *func_, int tightness_)
{
    if(eps >= 1)
        throw std::invalid_argument("eps >= 1");
    if(r_ <= 1)
        throw std::invalid_argument("r < 1");
    eps = eps_;
    r = r_;
    func = func_;
    iteration = 0;
    lowerBound = lowerBound_;
    upperBound = upperBound_;

    double lb[N], ub[N];
    for (size_t i = 0; i < N; ++i) {
        lb[i] = static_cast<double>(lowerBound_[i]);
        ub[i] = static_cast<double>(upperBound_[i]);
    }
    evolvent = ags::Evolvent(N, tightness_, lb, ub, ags::Simple);
}

template <class T, size_t N>
inline TAlgorithm<T, N>::~TAlgorithm()
{
}

template <class T, size_t N>
inline std::vector<T> TAlgorithm<T, N>::Solve(size_t maxInteration, bool isMinimize)
{
    T sign = isMinimize ? 1.0 : -1.0;
    double a = 0.0;
    double b = 1.0;
    double pointAData[N], pointBData[N];
    evolvent.GetImage(a, pointAData);
    evolvent.GetImage(b, pointBData);

    TPoint<T, N> pointA(pointAData), pointB(pointBData);

    T za = sign * func->ComputeFunction(pointA.toVector());
    T zb = sign * func->ComputeFunction(pointB.toVector());

    interval.initialize(a, b, za, zb);

    TPoint<T, N> bestPoint = pointA;
    T resZInternal = za;
    if (zb < resZInternal)
    {
        bestPoint = pointB;
        resZInternal = zb;
    }

    M = std::abs(zb - za) / pow((b - a), 1.0 / N);
    L = (M == 0) ? 1.0 : r * M;
    R = CalculateR(L, b, a, zb, za);
    interval.setIntervalR(0, R);

    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({R, 0});

    do
    {
        iteration++;
        
        indexInteravlWhithMaxR = pq.top().second;
        pq.pop();
        
        zi = interval.getZRight(indexInteravlWhithMaxR);
        zi1 = interval.getZLeft(indexInteravlWhithMaxR);
        xi = interval.getRight(indexInteravlWhithMaxR);
        xi1 = interval.getLeft(indexInteravlWhithMaxR);
        
        newX = ((xi1 + xi) / 2) - ((zi - zi1) / (2 * L));
        double newPointData[N];
        evolvent.GetImage(newX, newPointData);
        TPoint<T, N> newPoint (newPointData);
        T newZInternal = sign * func->ComputeFunction(newPoint.toVector());
        
        if (newZInternal < resZInternal)
        {
            bestPoint = newPoint;
            resZInternal = newZInternal;
        }

        size_t right_half_idx = interval.splitByIndex(indexInteravlWhithMaxR, newX, newZInternal);

        // Проверяем, выросла ли константа M
        double m1 = std::abs(newZInternal - zi1) / pow((newX - xi1), 1.0 / N);
        double m2 = std::abs(zi - newZInternal) / pow((xi - newX), 1.0 / N);
        double local_M = std::max(m1, m2);

        bool m_changed = false;
        if (local_M > M) 
        {
            M = local_M;
            L = (M == 0) ? 1.0 : r * M;
            m_changed = true;
        }

        if (m_changed) 
        {
            // т.к L изменилась пересчитываем все R
            pq = std::priority_queue<std::pair<double, size_t>>(); 
            for (size_t i = 0; i < interval.size(); i++)
            {
                double z_r = interval.getZRight(i);
                double z_l = interval.getZLeft(i);
                double x_r = interval.getRight(i);
                double x_l = interval.getLeft(i);
                
                double new_R = CalculateR(L, x_r, x_l, z_r, z_l);
                interval.setIntervalR(i, new_R);
                pq.push({new_R, i});
            }
        } 
        else 
        {
            // Пересчитываем R только для новых интервалов
            size_t idx1 = indexInteravlWhithMaxR; 
            size_t idx2 = right_half_idx;   
            
            double zl1 = interval.getZLeft(idx1);
            double zr1 = interval.getZRight(idx1);
            double xl1 = interval.getLeft(idx1);
            double xr1 = interval.getRight(idx1);
            double r1 = CalculateR(L, xr1, xl1, zr1, zl1);
            interval.setIntervalR(idx1, r1);
            pq.push({r1, idx1});
            
            double zl2 = interval.getZLeft(idx2);
            double zr2 = interval.getZRight(idx2);
            double xl2 = interval.getLeft(idx2);
            double xr2 = interval.getRight(idx2);
            double r2 = CalculateR(L, xr2, xl2, zr2, zl2);
            interval.setIntervalR(idx2, r2);
            pq.push({r2, idx2});
        }

    } 
    while ((interval.getLength(pq.top().second) > eps) && (iteration < maxInteration));

    size_t maxIndex = pq.top().second;
    double x0 = interval.getLeft(maxIndex);
    double x1 = interval.getRight(maxIndex);
    double xMid = (x0 + x1) / 2;
    double midPointData[N];
    evolvent.GetImage(xMid, midPointData);
    TPoint<T, N> midPoint;
    for (size_t i = 0; i < N; ++i)
        midPoint[i] = midPointData[i];
    
    T tmpZInternal = sign * func->ComputeFunction(midPoint.toVector());
    if (tmpZInternal < resZInternal)
    {
        bestPoint = midPoint;
        resZInternal = tmpZInternal;
    }
    
    T finalZ = isMinimize ? resZInternal : -resZInternal;
    std::vector<T> result;
    result.push_back(finalZ);
    for (size_t i = 0; i < N; ++i) 
        result.push_back(bestPoint[i]);
    result.push_back(static_cast<T>(iteration));
    return result;
};