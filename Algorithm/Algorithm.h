#include "Point.h"
#include "Interval.h"
#include "IOptProblem.hpp"
#include "Evolvent.hpp"
#include <queue>
#include <IGeneralOptProblem.hpp>

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
    IGeneralOptProblem* func;
    TInterval<double> interval;
    double eps;
    double r;
    double R;
    std::vector<double> M;
    std::vector<double> L;
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
    double CalculateR(const TInterval<double>& interval, const size_t index);
    double CalculateM(const TInterval<double>& interval, const size_t index);
    double CalculateNewX(const TInterval<double>& interval, const size_t index);
    void RebuildQueue(std::priority_queue<std::pair<double, size_t>>& pq, TInterval<double>& interval);
public:
    TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
              double r_, IGeneralOptProblem* func_, int tightness_);
    
    ~TAlgorithm();
    
    std::vector<T> Solve(size_t maxInteration, bool isMinimize = true); 
};


template <class T, size_t N>
inline double TAlgorithm<T, N>::CalculateR(const TInterval<double>& interval, const size_t index)
{
    double xLeft = interval.getLeft(index);
    double xRight = interval.getRight(index);
    double zLeft = interval.getZLeft(index);
    double zRight = interval.getZRight(index);
    double vLeft = interval.getVLeft(index);
    double vRight = interval.getVRight(index);
    double res = 0.0;
    if(vLeft == vRight)
    {
        res = pow((xRight - xLeft), 1.0 / N) + ((zRight - zLeft) * (zRight - zLeft)
        / (L[vLeft] * L[vLeft] * pow((xRight - xLeft), 1.0 / N))) - 2 * (zRight + zLeft) / L[vLeft];
    }
    else if(vLeft > vRight)
    {
        res = pow((xRight - xLeft), 1.0 / N) - 2 * zRight / L[vRight];
    }
    else
    {
        res = pow((xRight - xLeft), 1.0 / N) - 2 * zLeft / L[vLeft];
    }
   
    return res;
}

template <class T, size_t N>
inline double TAlgorithm<T, N>::CalculateM(const TInterval<double>& interval, const size_t index)
{
    double xLeft = interval.getLeft(index);
    double xRight = interval.getRight(index);
    double zLeft = interval.getZLeft(index);
    double zRight = interval.getZRight(index);
    double vLeft = interval.getVLeft(index);
    double vRight = interval.getVRight(index);

    if(vLeft == vRight)
    {
        double res = std::abs(zRight - zLeft) / pow((xRight - xLeft), 1.0 / N);
        return res;
    }
    return -1.0;
}

template <class T, size_t N>
inline double TAlgorithm<T, N>::CalculateNewX(const TInterval<double>& interval, const size_t index)
{
    double xLeft = interval.getLeft(index);
    double xRight = interval.getRight(index);
    double zLeft = interval.getZLeft(index);
    double zRight = interval.getZRight(index);
    double vLeft = interval.getVLeft(index);
    double vRight = interval.getVRight(index);

    double sgn = (zRight > zLeft) ? 1.0 : ((zRight < zLeft) ? -1.0 : 0.0);
    double newX = 0.0;
    if(vRight == vLeft)
    {
        newX = (xLeft + xRight) / 2.0 - sgn * pow(std::abs(zRight - zLeft) / L[vLeft], N) / 2.0;
    }
    else 
    {
        newX = (xLeft + xRight) / 2.0; 
    }
    return newX;
}

template <class T, size_t N>
inline void TAlgorithm<T, N>::RebuildQueue(std::priority_queue<std::pair<double, size_t>>& pq, TInterval<double> &interval)
{
    pq = std::priority_queue<std::pair<double, size_t>>();
    for (size_t i = 0; i < interval.size(); i++)
    {
        double new_R = CalculateR(interval, i);
        interval.setIntervalR(i, new_R);
        pq.push({new_R, i});
    }
}

template <class T, size_t N>
inline TAlgorithm<T, N>::TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
                                    double r_, IGeneralOptProblem* func_, int tightness_)
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
inline std::vector<T> TAlgorithm<T, N>::Solve(size_t maxIteration, bool isMinimize)
{
    T sign = isMinimize ? 1.0 : -1.0;
    // calculate first interval
    double ua = 0.0;
    double ub = 1.0;

    int evalIndex = func->GetConstraintsNumber();

    M.assign(evalIndex + 1, 0.0);
    L.assign(evalIndex + 1, 1.0);

    double pointAData[N], pointBData[N];
    evolvent.GetImage(ua, pointAData);
    evolvent.GetImage(ub, pointBData);
    int va = 0, vb = 0;
    TPoint<T, N> pointA(pointAData), pointB(pointBData);
    T za = sign * func->ComputePoint(std::vector<double>(pointAData, pointAData + N), va);
    T zb = sign * func->ComputePoint(std::vector<double>(pointBData, pointBData + N), vb);

    interval.initialize(ua, ub, za, zb, va, vb);

    TPoint<T, N> bestPoint = pointA;
    T resZInternal = za;
    if (vb == evalIndex) 
    {
        if (zb < resZInternal)
        {
            bestPoint = pointB;
            resZInternal = zb;
        }
    }

    M[va] = CalculateM(interval, 0);
    L[va] = (M[va] == 0) ? 1.0 : r * M[va];
    
    R = CalculateR(interval, 0);
    interval.setIntervalR(0, R);
    
    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({R, 0});
    do
    {
        iteration++;
        // get max R
        indexInteravlWhithMaxR = pq.top().second;
        pq.pop();

        // calculate a new point
        double newU = CalculateNewX(interval, indexInteravlWhithMaxR);
        int xEvalIndex = 0;
        // check value in point
        double newPointData[N];
        evolvent.GetImage(newU, newPointData);
        TPoint<T, N> newPoint(newPointData);
        double newZInternal = sign * func->ComputePoint(vector<double>(newPointData, newPointData + N), xEvalIndex);;

        if (xEvalIndex == evalIndex) 
        {
            if (newZInternal < resZInternal)
            {
                bestPoint = newPoint;
                resZInternal = newZInternal;
            }
        }
        // split the interval and get new two intervals

        size_t right_half_idx = interval.splitByIndex(indexInteravlWhithMaxR, newU, newZInternal, xEvalIndex);
        bool m_changed = false;
        // check has M changed
        // if yes, calculate for each new R 
        double m1 = CalculateM(interval, indexInteravlWhithMaxR);
        if((m1 > 0) && (M[xEvalIndex] < m1))
        {
            M[xEvalIndex] = m1;
            L[xEvalIndex] = (M[xEvalIndex] == 0) ? 1.0 : r * M[xEvalIndex];
            RebuildQueue(pq, interval);
            m_changed = true;
        }
        double m2 = CalculateM(interval, right_half_idx);
        if((m2 > 0) && (M[xEvalIndex] < m2))
        {
            M[xEvalIndex] = m2;
            L[xEvalIndex] = (M[xEvalIndex] == 0) ? 1.0 : r * M[xEvalIndex];
            RebuildQueue(pq, interval);
            m_changed = true;
        }

        // if no, calculate R for new intrefal
        if(!m_changed)
        {
            size_t idx1 = indexInteravlWhithMaxR;
            size_t idx2 = right_half_idx;

            
            double r1 = CalculateR(interval, idx1);
            interval.setIntervalR(idx1, r1);
            pq.push({r1, idx1});

            double r2 = CalculateR(interval, idx2);
            interval.setIntervalR(idx2, r2);
            pq.push({r2, idx2});
        }

    } while ((interval.getLength(pq.top().second) > eps) && (iteration < maxIteration));

    size_t maxIndex = pq.top().second;
    double u0 = interval.getLeft(maxIndex);
    double u1 = interval.getRight(maxIndex);
    double uMid = (u0 + u1) / 2;
    
    double midPointData[N];
    evolvent.GetImage(uMid, midPointData);
    
    TPoint<T, N> midPoint;
    for (size_t i = 0; i < N; ++i)
        midPoint[i] = midPointData[i];
    int tmp = 0;
    T tmpZInternal = sign * func->ComputePoint(vector<double>(midPointData, midPointData + N), tmp);
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
}