#include "Point.h"
#include "Interval.h"
#include "IOptProblem.hpp"
#include "Evolvent.hpp"
#include <queue>
#include <IGeneralOptProblem.hpp>
#include <climits>
#include "../../../Logger/Logger.h"

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
    std::vector<double> Z;
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
    double deltax = pow((xRight - xLeft), 1.0 / N);
    
    LOG_DEBUG("TAlgorithm::CalculateR - index={}, xLeft={}, xRight={}, zLeft={}, zRight={}, vLeft={}, vRight={}, deltax={}", 
              index, xLeft, xRight, zLeft, zRight, vLeft, vRight, deltax);
    
    if(vLeft == vRight)
    {
        res = deltax + ((zRight - zLeft) * (zRight - zLeft)
        / (L[vLeft] * L[vLeft] * deltax)) - 2 * ((zRight + zLeft) - 2 * Z[vLeft]) / L[vLeft];
        LOG_DEBUG("TAlgorithm::CalculateR - Case vLeft==vRight: res={}", res);
    }
    else if(vLeft > vRight)
    {
        res = 2 * deltax - 4 * (zRight - Z[vRight]) / L[vRight];
        LOG_DEBUG("TAlgorithm::CalculateR - Case vLeft>vRight: res={}", res);
    }
    else
    {
        res = 2 * deltax - 4 * (zLeft - Z[vLeft]) / L[vLeft];
        LOG_DEBUG("TAlgorithm::CalculateR - Case vLeft<vRight: res={}", res);
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

    LOG_DEBUG("TAlgorithm::CalculateM - index={}, xLeft={}, xRight={}, zLeft={}, zRight={}, vLeft={}, vRight={}", 
              index, xLeft, xRight, zLeft, zRight, vLeft, vRight);

    if(vLeft == vRight)
    {
        double res = std::abs(zRight - zLeft) / pow((xRight - xLeft), 1.0 / N);
        LOG_DEBUG("TAlgorithm::CalculateM - Case vLeft==vRight: res={}", res);
        return res;
    }
    
    LOG_DEBUG("TAlgorithm::CalculateM - Different vLeft and vRight, returning -1.0");
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

    LOG_DEBUG("TAlgorithm::CalculateNewX - index={}, xLeft={}, xRight={}, zLeft={}, zRight={}, vLeft={}, vRight={}", 
              index, xLeft, xRight, zLeft, zRight, vLeft, vRight);

    double dx = xRight - xLeft;
    if (dx <= 1e-15 || std::isnan(dx)) {
        LOG_DEBUG("TAlgorithm::CalculateNewX - Invalid dx: {}, returning midpoint", dx);
        return (xLeft + xRight) / 2.0;
    }

    double sgn = (zRight > zLeft) ? 1.0 : ((zRight < zLeft) ? -1.0 : 0.0);
    double newX = 0.0;
    if(vRight == vLeft && L[vLeft] > 0)
    {
        newX = (xLeft + xRight) / 2.0 - sgn * std::pow(std::abs(zRight - zLeft) / L[vLeft], static_cast<double>(N)) / 2.0;
        LOG_DEBUG("TAlgorithm::CalculateNewX - Calculated newX={}, sgn={}", newX, sgn);
        
        if (std::isnan(newX) || newX <= xLeft || newX >= xRight) {
            LOG_DEBUG("TAlgorithm::CalculateNewX - Invalid newX (NaN or out of bounds), using midpoint");
            newX = (xLeft + xRight) / 2.0;
        }
    }
    else 
    {
        LOG_DEBUG("TAlgorithm::CalculateNewX - vRight != vLeft or L[vLeft]<=0, using midpoint");
        newX = (xLeft + xRight) / 2.0; 
    }
    
    LOG_DEBUG("TAlgorithm::CalculateNewX - Final newX={}", newX);
    return newX;
}

template <class T, size_t N>
inline void TAlgorithm<T, N>::RebuildQueue(std::priority_queue<std::pair<double, size_t>>& pq, TInterval<double> &interval)
{
    LOG_DEBUG("TAlgorithm::RebuildQueue - Starting queue rebuild, current interval size={}", interval.size());
    pq = std::priority_queue<std::pair<double, size_t>>();
    for (size_t i = 0; i < interval.size(); i++)
    {
        double new_R = CalculateR(interval, i);
        interval.setIntervalR(i, new_R);
        pq.push({new_R, i});
    }
    LOG_DEBUG("TAlgorithm::RebuildQueue - Queue rebuild complete, new queue size={}", interval.size());
}

template <class T, size_t N>
inline TAlgorithm<T, N>::TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
                                    double r_, IGeneralOptProblem* func_, int tightness_)
{
    LOG_INFO("TAlgorithm::Constructor - Initializing with eps={}, r={}, tightness={}, N={}", 
             eps_, r_, tightness_, N);
    
    if(eps_ >= 1) {
        LOG_ERROR("TAlgorithm::Constructor - Invalid eps: {} (must be < 1)", eps_);
        throw std::invalid_argument("eps >= 1");
    }
    if(r_ <= 1) {
        LOG_ERROR("TAlgorithm::Constructor - Invalid r: {} (must be > 1)", r_);
        throw std::invalid_argument("r < 1");
    }
    
    eps = eps_;
    r = r_;
    func = func_;
    iteration = 0;
    lowerBound = lowerBound_;
    upperBound = upperBound_;

    LOG_DEBUG("TAlgorithm::Constructor - Lower bound: dim={}, Upper bound: dim={}", 
              N, N);

    double lb[N], ub[N];
    for (size_t i = 0; i < N; ++i) {
        lb[i] = static_cast<double>(lowerBound_[i]);
        ub[i] = static_cast<double>(upperBound_[i]);
    }
    evolvent = ags::Evolvent(N, tightness_, lb, ub, ags::Simple);
    
    LOG_INFO("TAlgorithm::Constructor - Successfully initialized");
}

template <class T, size_t N>
inline TAlgorithm<T, N>::~TAlgorithm()
{
    LOG_DEBUG("TAlgorithm::Destructor - Cleaning up");
}

template <class T, size_t N>
inline std::vector<T> TAlgorithm<T, N>::Solve(size_t maxIteration, bool isMinimize)
{
    LOG_INFO("TAlgorithm::Solve - Starting solve with maxIteration={}, isMinimize={}", maxIteration, isMinimize);
    
    T sign = isMinimize ? 1.0 : -1.0;
    // calculate first interval
    double ua = 0.0;
    double ub = 1.0;

    int evalIndex = func->GetConstraintsNumber();
    LOG_DEBUG("TAlgorithm::Solve - evalIndex (constraints number)={}", evalIndex);

    M.assign(evalIndex + 1, 0.0);
    L.assign(evalIndex + 1, 1.0);
    double MAX_double = std::numeric_limits<double>::max();
    Z.assign(evalIndex + 1, MAX_double);
    
    LOG_DEBUG("TAlgorithm::Solve - Initialized M, L, Z vectors with size={}", evalIndex + 1);
    
    double pointAData[N], pointBData[N];
    evolvent.GetImage(ua, pointAData);
    evolvent.GetImage(ub, pointBData);
    
    int va = 0, vb = 0;
    TPoint<T, N> pointA(pointAData), pointB(pointBData);
    T za = sign * func->ComputePoint(std::vector<double>(pointAData, pointAData + N), va);
    T zb = sign * func->ComputePoint(std::vector<double>(pointBData, pointBData + N), vb);

    LOG_DEBUG("TAlgorithm::Solve - Initial points: pointA z={}, va={}, pointB z={}, vb={}", za, va, zb, vb);

    interval.initialize(ua, ub, za, zb, va, vb);

    TPoint<T, N> bestPoint = pointA;
    T resZInternal = std::numeric_limits<T>::max();

    if (va == evalIndex) 
    {
        bestPoint = pointA;
        resZInternal = za;
        LOG_DEBUG("TAlgorithm::Solve - Best point found at A: z={}", za);
    }
    
    if (vb == evalIndex && zb < resZInternal) 
    {
        bestPoint = pointB;
        resZInternal = zb;
        LOG_DEBUG("TAlgorithm::Solve - Best point found at B: z={}", zb);
    }

    if(zb < Z[0])
    {
        Z[0] = zb;
        LOG_DEBUG("TAlgorithm::Solve - Updated Z[0] to {}", zb);
    }

    M[va] = CalculateM(interval, 0);
    L[va] = (M[va] == 0) ? 1.0 : r * M[va];
    LOG_DEBUG("TAlgorithm::Solve - Initial M[{}]={}, L[{}]={}", va, M[va], va, L[va]);
    
    R = CalculateR(interval, 0);
    interval.setIntervalR(0, R);
    LOG_DEBUG("TAlgorithm::Solve - Initial R={}", R);
    
    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({R, 0});
    
    do
    {
        iteration++;
        LOG_DEBUG("TAlgorithm::Solve - Iteration {}", iteration);
        
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
        double newZInternal = sign * func->ComputePoint(vector<double>(newPointData, newPointData + N), xEvalIndex);
        
        LOG_DEBUG("TAlgorithm::Solve - New point: u={}, z={}, evalIndex={}", newU, newZInternal, xEvalIndex);

        if (xEvalIndex == evalIndex) 
        {
            if (newZInternal < resZInternal)
            {
                bestPoint = newPoint;
                resZInternal = newZInternal;
                LOG_DEBUG("TAlgorithm::Solve - New best point found: z={}", resZInternal);
            }
        }
        if(xEvalIndex < Z.size() && newZInternal < Z[xEvalIndex])
        {
            Z[xEvalIndex] = newZInternal;
            LOG_DEBUG("TAlgorithm::Solve - Updated Z[{}] to {}", xEvalIndex, newZInternal);
        }
        // split the interval and get new two intervals

        size_t right_half_idx = interval.splitByIndex(indexInteravlWhithMaxR, newU, newZInternal, xEvalIndex);
        LOG_DEBUG("TAlgorithm::Solve - Split interval {} into {} and {}", indexInteravlWhithMaxR, indexInteravlWhithMaxR, right_half_idx);
        
        bool m_changed = false;
        // check has M changed
        // if yes, calculate for each new R 
        double m1 = CalculateM(interval, indexInteravlWhithMaxR);
        if(xEvalIndex < M.size() && (m1 > 0) && (M[xEvalIndex] < m1))
        {
            M[xEvalIndex] = m1;
            L[xEvalIndex] = (M[xEvalIndex] == 0) ? 1.0 : r * M[xEvalIndex];
            RebuildQueue(pq, interval);
            m_changed = true;
            LOG_DEBUG("TAlgorithm::Solve - M changed for xEvalIndex={}: new M={}, L={}", xEvalIndex, m1, L[xEvalIndex]);
        }
        double m2 = CalculateM(interval, right_half_idx);
        if(xEvalIndex < M.size() && (m2 > 0) && (M[xEvalIndex] < m2))
        {
            M[xEvalIndex] = m2;
            L[xEvalIndex] = (M[xEvalIndex] == 0) ? 1.0 : r * M[xEvalIndex];
            RebuildQueue(pq, interval);
            m_changed = true;
            LOG_DEBUG("TAlgorithm::Solve - M changed for xEvalIndex={}: new M={}, L={}", xEvalIndex, m2, L[xEvalIndex]);
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
            
            LOG_DEBUG("TAlgorithm::Solve - Added intervals to queue: R1={}, R2={}", r1, r2);
        }

        LOG_DEBUG("TAlgorithm::Solve - Queue size={}, Interval size={}", (int)interval.size(), (int)interval.size());

    } while ((interval.getLength(pq.top().second) > eps) && (iteration < maxIteration));

    LOG_INFO("TAlgorithm::Solve - Finished solve loop. Total iterations={}, final bestZ={}", iteration, resZInternal);

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
    
    LOG_DEBUG("TAlgorithm::Solve - Final midpoint z={}", tmpZInternal);
    
    if (tmpZInternal < resZInternal)
    {
        bestPoint = midPoint;
        resZInternal = tmpZInternal;
        LOG_DEBUG("TAlgorithm::Solve - Final point updated to midpoint: z={}", resZInternal);
    }

    T finalZ = isMinimize ? resZInternal : -resZInternal;
    std::vector<T> result;
    result.push_back(finalZ);
    for (size_t i = 0; i < N; ++i)
        result.push_back(bestPoint[i]);
    result.push_back(static_cast<T>(iteration));
    
    LOG_INFO("TAlgorithm::Solve - Solution: z={}, iterations={}", finalZ, iteration);
    return result;
}