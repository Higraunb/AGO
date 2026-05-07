#include "Point.h"
#include "Interval.h"
#include "IOptProblem.hpp"
#include "Evolvent.hpp"
#include "LocalMethod.hpp"
#include <queue>
#include <IGeneralOptProblem.hpp>
#include <climits>
#include "../../../Logger/Logger.h"
#include <algorithm>
#include <opencv2/ml.hpp>  


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

    std::vector<TPoint<T, N>> historyPoints;
    std::vector<double> historyValues;
    bool forestBuilt;
    size_t minTrialsBeforeForest;
    cv::Ptr<cv::ml::RTrees> randomForest; 

    double CalculateR(const TInterval<double>& interval, const size_t index);
    double CalculateM(const TInterval<double>& interval, const size_t index);
    double CalculateNewX(const TInterval<double>& interval, const size_t index);
    void RebuildQueue(std::priority_queue<std::pair<double, size_t>>& pq, TInterval<double>& interval);

    void BuildRandomForest();

public:
    TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
              double r_, IGeneralOptProblem* func_, int tightness_);
    
    ~TAlgorithm();
    
    std::vector<T> Solve0(size_t maxInteration, bool isMinimize = true);
    std::vector<T> Solve(size_t maxIteration, size_t p = 100, bool isMinimize = true);
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
    double deltax = pow((xRight - xLeft), 1.0 / N);
    LOG_DEBUG("TAlgorithm::CalculateM - index={}, xLeft={}, xRight={}, zLeft={}, zRight={}, vLeft={}, vRight={}", 
              index, xLeft, xRight, zLeft, zRight, vLeft, vRight);

    if(vLeft == vRight)
    {
        double res = std::abs(zRight - zLeft) / deltax;
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
    if(vRight == vLeft && M[vLeft] > 0)
    {
        newX = (xLeft + xRight) / 2.0 - sgn * std::pow(std::abs(zRight - zLeft) / M[vLeft], static_cast<double>(N)) / (2.0 * r);
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
inline void TAlgorithm<T, N>::BuildRandomForest()
{
    if (historyPoints.empty())
        return;

    LOG_DEBUG("TAlgorithm::BuildRandomForest - Building random forest with {} points", historyPoints.size());

    cv::Mat samples(static_cast<int>(historyPoints.size()), static_cast<int>(N), CV_32F);
    cv::Mat responses(static_cast<int>(historyPoints.size()), 1, CV_32F);

    for (size_t i = 0; i < historyPoints.size(); ++i) {
        for (size_t j = 0; j < N; ++j) {
            samples.at<float>(i, j) = static_cast<float>(historyPoints[i][j]);
        }
        responses.at<float>(i, 0) = static_cast<float>(historyValues[i]);
    }

    randomForest = cv::ml::RTrees::create();
    randomForest->setMaxDepth(6);
    randomForest->setMinSampleCount(10);
    randomForest->setRegressionAccuracy(0.01f);
    randomForest->setUseSurrogates(false);
    randomForest->setCVFolds(0);
    randomForest->setActiveVarCount(std::max(1, static_cast<int>(N / 2)));

    randomForest->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER, 30, 0.1));

    CV_Assert(samples.type() == CV_32F);
    CV_Assert(responses.type() == CV_32F);
    randomForest->train(samples, cv::ml::ROW_SAMPLE, responses);
    forestBuilt = true;

    LOG_DEBUG("TAlgorithm::BuildRandomForest - Random forest built successfully");
}

template<class T, size_t N>
inline TAlgorithm<T, N>::TAlgorithm(TPoint<T, N> lowerBound_, TPoint<T, N> upperBound_, double eps_,
                                    double r_, IGeneralOptProblem* func_, int tightness_)
    : func(func_),
      eps(eps_),
      r(r_),
      iteration(0),
      lowerBound(lowerBound_),
      upperBound(upperBound_),
      forestBuilt(false),
      minTrialsBeforeForest(100 * N)
{
    LOG_DEBUG("TAlgorithm::Constructor - Initializing with eps={}, r={}, tightness={}, N={}", 
              eps_, r_, tightness_, N);

    if (eps_ >= 1) {
        LOG_ERROR("TAlgorithm::Constructor - Invalid eps: {} (must be < 1)", eps_);
        throw std::invalid_argument("eps >= 1");
    }
    if (r_ <= 1) {
        LOG_ERROR("TAlgorithm::Constructor - Invalid r: {} (must be > 1)", r_);
        throw std::invalid_argument("r < 1");
    }

    LOG_DEBUG("TAlgorithm::Constructor - Lower bound: dim={}, Upper bound: dim={}", N, N);

    double lb[N], ub[N];
    for (size_t i = 0; i < N; ++i) {
        lb[i] = static_cast<double>(lowerBound_[i]);
        ub[i] = static_cast<double>(upperBound_[i]);
    }
    evolvent = ags::Evolvent(N, tightness_, lb, ub, ags::Simple);

    LOG_DEBUG("TAlgorithm::Constructor - Successfully initialized");
}

template <class T, size_t N>
inline TAlgorithm<T, N>::~TAlgorithm()
{
    LOG_DEBUG("TAlgorithm::Destructor - Cleaning up");
}

template <class T, size_t N>
inline std::vector<T> TAlgorithm<T, N>::Solve0(size_t maxIteration, bool isMinimize)
{
    LOG_DEBUG("TAlgorithm::Solve - Starting solve with maxIteration={}, isMinimize={}", maxIteration, isMinimize);
    
    T sign = isMinimize ? 1.0 : -1.0;
    // calculate first interval
    historyPoints.clear();
    historyValues.clear();
    forestBuilt = false;

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

    historyPoints.push_back(pointA);
    historyValues.push_back(static_cast<double>(za));
    historyPoints.push_back(pointB);
    historyValues.push_back(static_cast<double>(zb));

    TPoint<T, N> bestPoint = pointA;
    double resZInternal = std::min(za, zb);
    Z[0] = resZInternal;
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

        historyPoints.push_back(newPoint);
        historyValues.push_back(newZInternal);

        if (!forestBuilt && historyPoints.size() >= minTrialsBeforeForest) {
            BuildRandomForest();
            //break;
        }

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

    LOG_DEBUG("TAlgorithm::Solve - Finished solve loop. Total iterations={}, final bestZ={}", iteration, resZInternal);

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

    historyPoints.push_back(midPoint);
    historyValues.push_back(static_cast<double>(tmpZInternal));
    
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
    
    LOG_DEBUG("TAlgorithm::Solve - Solution: z={}, iterations={}", finalZ, iteration);
    return result;
}

template <class T, size_t N>
inline std::vector<T> TAlgorithm<T, N>::Solve(size_t maxIteration, size_t p, bool isMinimize)
{
    LOG_DEBUG("TAlgorithm::Solve - Starting solve with maxIteration={}, p={}, isMinimize={}", maxIteration, p, isMinimize);
    
    T sign = isMinimize ? 1.0 : -1.0;
    historyPoints.clear();
    historyValues.clear();
    forestBuilt = false;

    double ua = 0.0;
    double ub = 1.0;

    int evalIndex = func->GetConstraintsNumber();
    M.assign(evalIndex + 1, 0.0);
    L.assign(evalIndex + 1, 1.0);
    double MAX_double = std::numeric_limits<double>::max();
    Z.assign(evalIndex + 1, MAX_double);
    
    double pointAData[N], pointBData[N];
    evolvent.GetImage(ua, pointAData);
    evolvent.GetImage(ub, pointBData);
    
    int va = 0, vb = 0;
    TPoint<T, N> pointA(pointAData), pointB(pointBData);
    T za = sign * func->ComputePoint(std::vector<double>(pointAData, pointAData + N), va);
    T zb = sign * func->ComputePoint(std::vector<double>(pointBData, pointBData + N), vb);

    interval.initialize(ua, ub, za, zb, va, vb);

    historyPoints.push_back(pointA);
    historyValues.push_back(static_cast<double>(za));
    historyPoints.push_back(pointB);
    historyValues.push_back(static_cast<double>(zb));

    TPoint<T, N> bestPoint = pointA;
    double resZInternal = std::min(za, zb);
    Z[0] = resZInternal;
    M[va] = CalculateM(interval, 0);
    L[va] = (M[va] == 0) ? 1.0 : r * M[va];
    
    R = CalculateR(interval, 0);
    interval.setIntervalR(0, R);
    
    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({R, 0});

    size_t lastTrainedSize = 0;
    
    do
    {
        iteration++;
        
        indexInteravlWhithMaxR = pq.top().second;
        pq.pop();

        double newU = CalculateNewX(interval, indexInteravlWhithMaxR);
        int xEvalIndex = 0;
        double newPointData[N];
        evolvent.GetImage(newU, newPointData);
        TPoint<T, N> newPoint(newPointData);
        double newZInternal = sign * func->ComputePoint(std::vector<double>(newPointData, newPointData + N), xEvalIndex);

        historyPoints.push_back(newPoint);
        historyValues.push_back(newZInternal);

        if (!forestBuilt && historyPoints.size() >= minTrialsBeforeForest) {
            BuildRandomForest();
        }

        if (xEvalIndex == evalIndex && newZInternal < resZInternal) {
            bestPoint = newPoint;
            resZInternal = newZInternal;
        }
        
        if(xEvalIndex < Z.size() && newZInternal < Z[xEvalIndex]) {
            Z[xEvalIndex] = newZInternal;
        }

        size_t right_half_idx = interval.splitByIndex(indexInteravlWhithMaxR, newU, newZInternal, xEvalIndex);
        
        bool m_changed = false;
        double m1 = CalculateM(interval, indexInteravlWhithMaxR);
        if(xEvalIndex < M.size() && (m1 > 0) && (M[xEvalIndex] < m1)) {
            M[xEvalIndex] = m1;
            L[xEvalIndex] = (M[xEvalIndex] == 0) ? 1.0 : r * M[xEvalIndex];
            RebuildQueue(pq, interval);
            m_changed = true;
        }
        double m2 = CalculateM(interval, right_half_idx);
        if(xEvalIndex < M.size() && (m2 > 0) && (M[xEvalIndex] < m2)) {
            M[xEvalIndex] = m2;
            L[xEvalIndex] = (M[xEvalIndex] == 0) ? 1.0 : r * M[xEvalIndex];
            if(!m_changed) RebuildQueue(pq, interval);
            m_changed = true;
        }

        if(!m_changed) {
            double r1 = CalculateR(interval, indexInteravlWhithMaxR);
            interval.setIntervalR(indexInteravlWhithMaxR, r1);
            pq.push({r1, indexInteravlWhithMaxR});

            double r2 = CalculateR(interval, right_half_idx);
            interval.setIntervalR(right_half_idx, r2);
            pq.push({r2, right_half_idx});
        }

        // --- ВНЕДРЕНИЕ: СЕТКА, RF И ЛОКАЛЬНЫЙ ПОИСК ---
        static size_t current_p = p; // Используем динамический шаг

        if (historyPoints.size() - lastTrainedSize >= current_p) {
            BuildRandomForest();
            
            static size_t lastCheckedSize = 0; 
            size_t startIdx = lastCheckedSize;
            size_t endIdx = historyPoints.size();
            
            lastTrainedSize = historyPoints.size();
            lastCheckedSize = historyPoints.size();

            // 1. Динамически увеличиваем шаг p, чтобы реже переобучать лес на глубоких этапах
            current_p = static_cast<size_t>(current_p * 1.5); 

            // 2. Вычисляем шаги (дельты) для "соседей по сетке"
            const int GRID_RES = 15;
            std::vector<double> delta(N);
            for (size_t d = 0; d < N; ++d) {
                delta[d] = (static_cast<double>(upperBound[d]) - static_cast<double>(lowerBound[d])) / (GRID_RES - 1.0);
            }

            // 3. ФИЛЬТРАЦИЯ "БЕСПЕРСПЕКТИВНЫХ" ТОЧЕК
            // Нет смысла искать локальный экстремум там, где значение функции заведомо хуже рекорда
            std::vector<std::pair<double, size_t>> recentPoints;
            for (size_t i = startIdx; i < endIdx; ++i) {
                recentPoints.push_back({historyValues[i], i});
            }

            // Сортируем: лучшие значения оказываются в начале
            if (isMinimize) {
                std::sort(recentPoints.begin(), recentPoints.end(), 
                    [](const auto& a, const auto& b) { return a.first < b.first; });
            } else {
                std::sort(recentPoints.begin(), recentPoints.end(), 
                    [](const auto& a, const auto& b) { return a.first > b.first; });
            }

            // Берем максимум 3 самые перспективные точки из новых
            size_t topK = std::min<size_t>(3, recentPoints.size());
            std::vector<size_t> promisingIndices;
            for (size_t i = 0; i < topK; ++i) {
                promisingIndices.push_back(recentPoints[i].second);
            }

            bool foundAnyCandidate = false;
            TPoint<T, N> overallBestLocalPoint;
            double overallBestLocalZ = isMinimize ? std::numeric_limits<double>::max() : -std::numeric_limits<double>::max();

            if (!promisingIndices.empty()) {
                size_t numNeighbors = 1ULL << N; // 2^n соседей
                
                // 4. ПАКЕТНОЕ СОЗДАНИЕ ВЫБОРКИ И ПРЕДСКАЗАНИЕ (BATCH PREDICTION)
                cv::Mat neighborSamples(promisingIndices.size() * numNeighbors, (int)N, CV_32F);
                
                size_t rowIdx = 0;
                for (size_t i : promisingIndices) {
                    const auto& pt = historyPoints[i];
                    for (size_t j = 0; j < numNeighbors; ++j) {
                        for (size_t d = 0; d < N; ++d) {
                            // Битовая маска для генерации 2^n направлений
                            double shift = ((j >> d) & 1) ? delta[d] : -delta[d];
                            double neighborVal = static_cast<double>(pt[d]) + shift;
                            neighborVal = std::max(static_cast<double>(lowerBound[d]), 
                                          std::min(static_cast<double>(upperBound[d]), neighborVal));
                            neighborSamples.at<float>(rowIdx, d) = static_cast<float>(neighborVal);
                        }
                        rowIdx++;
                    }
                }

                // Один вызов предсказания для всех соседей всех перспективных точек
                cv::Mat predictions;
                randomForest->predict(neighborSamples, predictions);

                // 5. ПОИСК ЛОКАЛЬНЫХ ЭКСТРЕМУМОВ И ЗАПУСК ХУКА-ДЖИВСА
                for (size_t idx = 0; idx < promisingIndices.size(); ++idx) {
                    size_t ptIndex = promisingIndices[idx];
                    double ptVal = historyValues[ptIndex];
                    bool isLocalExtremum = true;
                    
                    size_t baseRow = idx * numNeighbors;
                    for (size_t j = 0; j < numNeighbors; ++j) {
                        double predZ = predictions.at<float>(baseRow + j, 0);
                        if (isMinimize ? (ptVal >= predZ) : (ptVal <= predZ)) { 
                            isLocalExtremum = false; 
                            break; 
                        }
                    }

                    if (isLocalExtremum) {
                        LocalMethod<T, N> lm(func, historyPoints[ptIndex], lowerBound, upperBound, isMinimize);
                        lm.SetMaxTrial(40); // Снижаем лимит для скорости (раньше было 70)
                        auto localRes = lm.StartOptimization();
                        
                        double internalLocalZ = isMinimize ? static_cast<double>(localRes.second) : -static_cast<double>(localRes.second);

                        if (!foundAnyCandidate || 
                           (isMinimize ? (internalLocalZ < overallBestLocalZ) : (internalLocalZ > overallBestLocalZ))) {
                            overallBestLocalZ = internalLocalZ;
                            overallBestLocalPoint = localRes.first;
                            foundAnyCandidate = true;
                        }
                    }
                }
            }

            // 6. Интеграция лучшей локальной точки обратно в АГП
            if (foundAnyCandidate) {
                int xEvalLocalIndex = 0;
                std::vector<double> bestLocalDouble(N);
                for(size_t d = 0; d < N; d++) bestLocalDouble[d] = static_cast<double>(overallBestLocalPoint[d]);
                
                double trueZLocalInternal = sign * func->ComputePoint(bestLocalDouble, xEvalLocalIndex);

                if (isMinimize ? (trueZLocalInternal < resZInternal) : (trueZLocalInternal > resZInternal)) {
                    if (xEvalLocalIndex == evalIndex) {
                        resZInternal = trueZLocalInternal;
                        bestPoint = overallBestLocalPoint;
                    }
                }

                if (xEvalLocalIndex < Z.size() && trueZLocalInternal < Z[xEvalLocalIndex]) {
                    Z[xEvalLocalIndex] = trueZLocalInternal;
                }

                historyPoints.push_back(overallBestLocalPoint);
                historyValues.push_back(trueZLocalInternal);

                double u_preimages[ags::noninjectiveMaxPreimages]; 
                int preimCount = evolvent.GetAllPreimages(bestLocalDouble.data(), u_preimages);
                
                if (preimCount > 0) {
                    double u_new = u_preimages[0]; 
                    size_t targetIntervalIndex = std::numeric_limits<size_t>::max();
                    for (size_t i = 0; i < interval.size(); ++i) {
                        if (u_new > interval.getLeft(i) && u_new < interval.getRight(i)) {
                            targetIntervalIndex = i;
                            break;
                        }
                    }

                    if (targetIntervalIndex != std::numeric_limits<size_t>::max()) {
                        size_t right_half_idx = interval.splitByIndex(targetIntervalIndex, u_new, trueZLocalInternal, xEvalLocalIndex);
                        
                        bool m_changed = false;
                        double m1 = CalculateM(interval, targetIntervalIndex);
                        if(xEvalLocalIndex < M.size() && (m1 > 0) && (M[xEvalLocalIndex] < m1)) {
                            M[xEvalLocalIndex] = m1;
                            L[xEvalLocalIndex] = (M[xEvalLocalIndex] == 0) ? 1.0 : r * M[xEvalLocalIndex];
                            RebuildQueue(pq, interval);
                            m_changed = true;
                        }
                        double m2 = CalculateM(interval, right_half_idx);
                        if(xEvalLocalIndex < M.size() && (m2 > 0) && (M[xEvalLocalIndex] < m2)) {
                            M[xEvalLocalIndex] = m2;
                            L[xEvalLocalIndex] = (M[xEvalLocalIndex] == 0) ? 1.0 : r * M[xEvalLocalIndex];
                            if(!m_changed) RebuildQueue(pq, interval);
                            m_changed = true;
                        }

                        if(!m_changed) {
                            double r1 = CalculateR(interval, targetIntervalIndex);
                            interval.setIntervalR(targetIntervalIndex, r1);
                            pq.push({r1, targetIntervalIndex});

                            double r2 = CalculateR(interval, right_half_idx);
                            interval.setIntervalR(right_half_idx, r2);
                            pq.push({r2, right_half_idx});
                        }
                    }
                }
            }
        }

    } while ((interval.getLength(pq.top().second) > eps) && (iteration < maxIteration));

    T finalZ = isMinimize ? resZInternal : -resZInternal;
    std::vector<T> result;
    result.push_back(finalZ);
    for (size_t i = 0; i < N; ++i) result.push_back(bestPoint[i]);
    result.push_back(static_cast<T>(iteration));
    
    return result;
}
