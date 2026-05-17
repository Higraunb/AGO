#pragma once

#include "Point.h"
#include "Interval.h"
#include "IOptProblem.hpp"
#include "Evolvent.hpp"
#include "LocalMethod.hpp"
#include <queue>
#include <vector>
#include <climits>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <opencv2/ml.hpp>
#include "../../../Logger/Logger.h"
#include <IGeneralOptProblem.hpp>

template<size_t N>
class Algorithm
{
private:
    IGeneralOptProblem* func;
    Intervals interval; 
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
    TPoint<N> lowerBound;
    TPoint<N> upperBound;
    size_t indexInteravlWhithMaxR;
    size_t indexInteravlWhithMinR;
    size_t iteration;

    std::vector<TPoint<N>> historyPoints; 
    std::vector<double> historyValues;
    bool forestBuilt;
    size_t minTrialsBeforeForest;
    cv::Ptr<cv::ml::RTrees> randomForest; 

    double CalculateR(const Intervals& interval, const size_t index);
    double CalculateM(const Intervals& interval, const size_t index);
    double CalculateNewX(const Intervals& interval, const size_t index);
    void RebuildQueue(std::priority_queue<std::pair<double, size_t>>& pq, Intervals& interval);

    void BuildRandomForest();

    inline double FastPowInv(double dx) 
    {
        if constexpr (N == 1) return dx;
        else if constexpr (N == 2) return std::sqrt(dx);
        else if constexpr (N == 3) return std::cbrt(dx);
        else if constexpr (N == 4) return std::sqrt(std::sqrt(dx));
        else return std::pow(dx, 1.0 / static_cast<double>(N));
    }

    inline double FastPowInt(double base) 
    {
        double res = 1.0;
        for (size_t i = 0; i < N; ++i) 
        {
            res *= base;
        }
        return res;
    }

public:
    Algorithm(const TPoint<N>& lowerBound_, const TPoint<N>& upperBound_, double eps_,
              double r_, IGeneralOptProblem* func_, int tightness_);
    
    ~Algorithm();
    
    std::vector<double> Solve0(size_t maxIteration, bool isMinimize = true);
    std::vector<double> Solve(size_t maxIteration, size_t p = 100, bool isMinimize = true);
};

template <size_t N>
inline double Algorithm<N>::CalculateR(const Intervals& interval, const size_t index)
{
    double xLeft = interval.getLeft(index);
    double xRight = interval.getRight(index);
    double zLeft = interval.getZLeft(index);
    double zRight = interval.getZRight(index);
    double vLeft = interval.getVLeft(index);
    double vRight = interval.getVRight(index);
    double res = 0.0;
    double deltax = FastPowInv(xRight - xLeft);
    
    if (vLeft == vRight)
    {
        res = deltax + ((zRight - zLeft) * (zRight - zLeft)
        / (L[vLeft] * L[vLeft] * deltax)) - 2 * ((zRight + zLeft) - 2 * Z[vLeft]) / L[vLeft];
    }
    else if (vLeft > vRight)
    {
        res = 2 * deltax - 4 * (zRight - Z[vRight]) / L[vRight];
    }
    else
    {
        res = 2 * deltax - 4 * (zLeft - Z[vLeft]) / L[vLeft];
    }

    return res;
}

template <size_t N>
inline double Algorithm<N>::CalculateM(const Intervals& interval, const size_t index)
{
    double xLeft = interval.getLeft(index);
    double xRight = interval.getRight(index);
    double zLeft = interval.getZLeft(index);
    double zRight = interval.getZRight(index);
    double vLeft = interval.getVLeft(index);
    double vRight = interval.getVRight(index);
    double deltax = FastPowInv(xRight - xLeft);

    if (vLeft == vRight)
    {
        double res = std::abs(zRight - zLeft) / deltax;
        return res;
    }
    
    return -1.0;
}

template <size_t N>
inline double Algorithm<N>::CalculateNewX(const Intervals& interval, const size_t index)
{
    double xLeft = interval.getLeft(index);
    double xRight = interval.getRight(index);
    double zLeft = interval.getZLeft(index);
    double zRight = interval.getZRight(index);
    double vLeft = interval.getVLeft(index);
    double vRight = interval.getVRight(index);

    double dx = xRight - xLeft;
    if (dx <= 1e-15 || std::isnan(dx)) 
    {
        return (xLeft + xRight) / 2.0;
    }

    double sgn = (zRight > zLeft) ? 1.0 : ((zRight < zLeft) ? -1.0 : 0.0);
    double newX = 0.0;
    if (vRight == vLeft && M[vLeft] > 0)
    {
        double tmp = std::abs(zRight - zLeft) / M[vLeft];
        newX = (xLeft + xRight) / 2.0 - sgn * FastPowInt(tmp) / (2.0 * r);
       
        if (std::isnan(newX) || newX <= xLeft || newX >= xRight) 
        {
            newX = (xLeft + xRight) / 2.0;
        }
    }
    else 
    {
        newX = (xLeft + xRight) / 2.0; 
    }
 
    return newX;
}

template <size_t N>
inline void Algorithm<N>::RebuildQueue(std::priority_queue<std::pair<double, size_t>>& pq, Intervals &interval)
{
    std::vector<std::pair<double, size_t>> container;
    container.reserve(interval.size());
    for (size_t i = 0; i < interval.size(); i++)
    {
        double new_R = CalculateR(interval, i);
        interval.setIntervalR(i, new_R);
        container.push_back({new_R, i});
    }
    
    std::priority_queue<std::pair<double, size_t>> new_pq(
        std::less<std::pair<double, size_t>>(), 
        std::move(container)
    );
    pq.swap(new_pq);
}

template <size_t N>
inline void Algorithm<N>::BuildRandomForest()
{
    if (historyPoints.empty())
    {
        return;
    }

    const size_t MAX_SAMPLES = 10000;
    size_t start_idx = (historyPoints.size() > MAX_SAMPLES) ? historyPoints.size() - MAX_SAMPLES : 0;
    size_t sample_count = historyPoints.size() - start_idx;

    cv::Mat samples(static_cast<int>(sample_count), static_cast<int>(N), CV_32F);
    cv::Mat responses(static_cast<int>(sample_count), 1, CV_32F);

    for (size_t i = 0; i < sample_count; ++i) 
    {
        size_t actual_idx = start_idx + i;
        for (size_t j = 0; j < N; ++j) 
        {
            samples.at<float>(i, j) = static_cast<float>(historyPoints[actual_idx][j]);
        }
        responses.at<float>(i, 0) = static_cast<float>(historyValues[actual_idx]);
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
}

template <size_t N>
inline Algorithm<N>::Algorithm(const TPoint<N>& lowerBound_, const TPoint<N>& upperBound_, double eps_,
                            double r_, IGeneralOptProblem* func_, int tightness_)
    : func(func_),
      eps(eps_),
      r(r_),
      iteration(0),
      lowerBound(lowerBound_),
      upperBound(upperBound_),
      forestBuilt(false)
{
    minTrialsBeforeForest = 100 * N;
    
    if (eps_ >= 1) throw std::invalid_argument("eps >= 1");
    if (r_ <= 1) throw std::invalid_argument("r < 1");

    evolvent = ags::Evolvent(N, tightness_, const_cast<double*>(&lowerBound[0]), const_cast<double*>(&upperBound[0]), ags::Simple);
}

template <size_t N>
inline Algorithm<N>::~Algorithm()
{
}

template <size_t N>
inline std::vector<double> Algorithm<N>::Solve0(size_t maxIteration, bool isMinimize)
{
    double sign = isMinimize ? 1.0 : -1.0;
    
    historyPoints.clear();
    historyPoints.reserve(maxIteration + 100);
    historyValues.reserve(maxIteration + 100);
    historyValues.clear();
    forestBuilt = false;

    double ua = 0.0;
    double ub = 1.0;

    int evalIndex = func->GetConstraintsNumber();

    M.assign(evalIndex + 1, 0.0);
    L.assign(evalIndex + 1, 1.0);
    double MAX_double = std::numeric_limits<double>::max();
    Z.assign(evalIndex + 1, MAX_double);
    
    double pointAData[N];
    double pointBData[N];
    evolvent.GetImage(ua, pointAData);
    evolvent.GetImage(ub, pointBData);
    
    int va = 0, vb = 0;
    TPoint<N> pointA;
    TPoint<N> pointB;
    std::vector<double> vecA(N);
    std::vector<double> vecB(N);

    for (size_t i = 0; i < N; ++i) 
    {
        pointA[i] = pointAData[i];
        pointB[i] = pointBData[i];
        vecA[i] = pointAData[i];
        vecB[i] = pointBData[i];
    }
    
    double za = sign * func->ComputePoint(vecA, va); 
    double zb = sign * func->ComputePoint(vecB, vb);

    interval.initialize(ua, ub, za, zb, va, vb);

    historyPoints.push_back(pointA);
    historyValues.push_back(static_cast<double>(za));
    historyPoints.push_back(pointB);
    historyValues.push_back(static_cast<double>(zb));

    TPoint<N> bestPoint = pointA;
    double resZInternal = std::min(za, zb);
    Z[0] = resZInternal;
    M[va] = CalculateM(interval, 0);
    L[va] = (M[va] == 0) ? 1.0 : r * M[va];
   
    R = CalculateR(interval, 0);
    interval.setIntervalR(0, R);
    
    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({R, 0});
    double resPoint = func->GetOptimumValue();
    do
    {
        iteration++;
       
        indexInteravlWhithMaxR = pq.top().second;
        pq.pop();

        double newU = CalculateNewX(interval, indexInteravlWhithMaxR);
        int xEvalIndex = 0;
        
        double newPointData[N];
        evolvent.GetImage(newU, newPointData);
        TPoint<N> newPoint;
        std::vector<double> newPointVec(N);

        for (size_t i = 0; i < N; ++i) 
        {
            newPoint[i] = newPointData[i];
            newPointVec[i] = newPointData[i];
        }

        double newZInternal = sign * func->ComputePoint(newPointVec, xEvalIndex);

        historyPoints.push_back(newPoint);
        historyValues.push_back(newZInternal);

        if (!forestBuilt && historyPoints.size() >= minTrialsBeforeForest) 
        {
            BuildRandomForest();
        }

        if (xEvalIndex == evalIndex) 
        {
            if (newZInternal < resZInternal)
            {
                bestPoint = newPoint;
                resZInternal = newZInternal;
            }
        }
        
        if (xEvalIndex < Z.size() && newZInternal < Z[xEvalIndex])
        {
            Z[xEvalIndex] = newZInternal;
        }

        size_t right_half_idx = interval.splitByIndex(indexInteravlWhithMaxR, newU, newZInternal, xEvalIndex);
        
        bool m_changed = false;
        
        double m1 = CalculateM(interval, indexInteravlWhithMaxR);
        if (xEvalIndex < M.size() && (m1 > 0) && (M[xEvalIndex] < m1))
        {
            M[xEvalIndex] = m1;
            L[xEvalIndex] = (M[xEvalIndex] == 0) ? 1.0 : r * M[xEvalIndex];
            RebuildQueue(pq, interval);
            m_changed = true;
        }
        
        double m2 = CalculateM(interval, right_half_idx);
        if (xEvalIndex < M.size() && (m2 > 0) && (M[xEvalIndex] < m2))
        {
            M[xEvalIndex] = m2;
            L[xEvalIndex] = (M[xEvalIndex] == 0) ? 1.0 : r * M[xEvalIndex];
            RebuildQueue(pq, interval);
            m_changed = true;
        }

        if (!m_changed)
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

    } while ((std::abs(resZInternal - resPoint) > eps) && (iteration < maxIteration));


    size_t maxIndex = pq.top().second;
    double u0 = interval.getLeft(maxIndex);
    double u1 = interval.getRight(maxIndex);
    double uMid = (u0 + u1) / 2.0;
    
    double midPointData[N];
    evolvent.GetImage(uMid, midPointData);
    
    TPoint<N> midPoint;
    std::vector<double> midPointVec(N);

    for (size_t i = 0; i < N; ++i) 
    {
        midPoint[i] = midPointData[i];
        midPointVec[i] = midPointData[i];
    }

    int tmp = 0;
    double tmpZInternal = sign * func->ComputePoint(midPointVec, tmp);

    historyPoints.push_back(midPoint);
    historyValues.push_back(static_cast<double>(tmpZInternal));
    
    if (tmpZInternal < resZInternal)
    {
        bestPoint = midPoint;
        resZInternal = tmpZInternal;
    }

    double finalZ = isMinimize ? resZInternal : -resZInternal;

    std::vector<double> result;
    result.push_back(finalZ);
    for (size_t i = 0; i < N; ++i)
    {
        result.push_back(bestPoint[i]);
    }
    result.push_back(static_cast<double>(iteration));
    
    return result;
}

template <size_t N>
inline std::vector<double> Algorithm<N>::Solve(size_t maxIteration, size_t p, bool isMinimize)
{
    std::vector<double> bestLocalDouble(N);
    std::vector<double> delta(N);

    double sign = isMinimize ? 1.0 : -1.0;
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
    
    double pointAData[N];
    double pointBData[N];
    evolvent.GetImage(ua, pointAData);
    evolvent.GetImage(ub, pointBData);
    
    int va = 0, vb = 0;
    TPoint<N> pointA;
    TPoint<N> pointB;
    std::vector<double> vecA(N);
    std::vector<double> vecB(N);

    for (size_t i = 0; i < N; ++i) 
    {
        pointA[i] = pointAData[i];
        pointB[i] = pointBData[i];
        vecA[i] = pointAData[i];
        vecB[i] = pointBData[i];
    }

    double za = sign * func->ComputePoint(vecA, va);
    double zb = sign * func->ComputePoint(vecB, vb);

    interval.initialize(ua, ub, za, zb, va, vb);

    historyPoints.push_back(pointA);
    historyValues.push_back(static_cast<double>(za));
    historyPoints.push_back(pointB);
    historyValues.push_back(static_cast<double>(zb));

    TPoint<N> bestPoint = pointA;
    double resZInternal = std::min(za, zb);
    Z[0] = resZInternal;
    M[va] = CalculateM(interval, 0);
    L[va] = (M[va] == 0) ? 1.0 : r * M[va];
    
    R = CalculateR(interval, 0);
    interval.setIntervalR(0, R);
    
    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({R, 0});

    const size_t MAX_TOP_K = 3;
    const size_t numNeighbors = 1ULL << N; 
    cv::Mat neighborSamplesBuffer(static_cast<int>(MAX_TOP_K * numNeighbors), static_cast<int>(N), CV_32F);
    cv::Mat predictions;

    size_t lastTrainedSize = 0;
    size_t lastCheckedSize = 0; 
    size_t current_p = p;
    double resPoint = func->GetOptimumValue();
    do
    {
        iteration++;
        
        indexInteravlWhithMaxR = pq.top().second;
        pq.pop();

        double newU = CalculateNewX(interval, indexInteravlWhithMaxR);
        int xEvalIndex = 0;
        
        double newPointData[N];
        evolvent.GetImage(newU, newPointData);
        TPoint<N> newPoint;
        std::vector<double> newPointVec(N);

        for (size_t i = 0; i < N; ++i) 
        {
            newPoint[i] = newPointData[i];
            newPointVec[i] = newPointData[i];
        }
        
        double newZInternal = sign * func->ComputePoint(newPointVec, xEvalIndex);

        historyPoints.push_back(newPoint);
        historyValues.push_back(newZInternal);

        if (!forestBuilt && historyPoints.size() >= minTrialsBeforeForest) 
        {
            BuildRandomForest();
        }

        if (xEvalIndex == evalIndex && newZInternal < resZInternal) 
        {
            bestPoint = newPoint;
            resZInternal = newZInternal;
        }
        
        if (xEvalIndex < Z.size() && newZInternal < Z[xEvalIndex])
        {
            Z[xEvalIndex] = newZInternal;
        }

        size_t right_half_idx = interval.splitByIndex(indexInteravlWhithMaxR, newU, newZInternal, xEvalIndex);
        
        bool m_changed = false;
        
        double m1 = CalculateM(interval, indexInteravlWhithMaxR);
        if (xEvalIndex < M.size() && (m1 > 0) && (M[xEvalIndex] < m1)) 
        {
            M[xEvalIndex] = m1;
            L[xEvalIndex] = (M[xEvalIndex] == 0) ? 1.0 : r * M[xEvalIndex];
            m_changed = true;
        }
        
        double m2 = CalculateM(interval, right_half_idx);
        if (xEvalIndex < M.size() && (m2 > 0) && (M[xEvalIndex] < m2)) 
        {
            M[xEvalIndex] = m2;
            L[xEvalIndex] = (M[xEvalIndex] == 0) ? 1.0 : r * M[xEvalIndex];
            m_changed = true;
        }

        if (m_changed)
        {
            RebuildQueue(pq, interval);
        }
        else 
        {
            double r1 = CalculateR(interval, indexInteravlWhithMaxR);
            interval.setIntervalR(indexInteravlWhithMaxR, r1);
            pq.push({r1, indexInteravlWhithMaxR});

            double r2 = CalculateR(interval, right_half_idx);
            interval.setIntervalR(right_half_idx, r2);
            pq.push({r2, right_half_idx});
        }

        if (historyPoints.size() - lastTrainedSize >= current_p) 
        {
            BuildRandomForest();
            
            size_t startIdx = lastCheckedSize;
            size_t endIdx = historyPoints.size();
            
            lastTrainedSize = historyPoints.size();
            lastCheckedSize = historyPoints.size();

            current_p = static_cast<size_t>(current_p * 1.5); 

            const int GRID_RES = 15;
            
            for (size_t d = 0; d < N; ++d) 
            {
                delta[d] = (upperBound[d] - lowerBound[d]) / (GRID_RES - 1.0);
            }

            std::vector<std::pair<double, size_t>> recentPoints;
            for (size_t i = startIdx; i < endIdx; ++i) 
            {
                recentPoints.push_back({historyValues[i], i});
            }

            if (isMinimize) 
            {
                std::sort(recentPoints.begin(), recentPoints.end(), 
                    [](const auto& a, const auto& b) { return a.first < b.first; });
            } 
            else 
            {
                std::sort(recentPoints.begin(), recentPoints.end(), 
                    [](const auto& a, const auto& b) { return a.first > b.first; });
            }

            size_t topK = std::min<size_t>(3, recentPoints.size());
            std::vector<size_t> promisingIndices;
            for (size_t i = 0; i < topK; ++i) 
            {
                promisingIndices.push_back(recentPoints[i].second);
            }

            bool foundAnyCandidate = false;
            TPoint<N> overallBestLocalPoint;
            double overallBestLocalZ = isMinimize ? std::numeric_limits<double>::max() : -std::numeric_limits<double>::max();

            if (!promisingIndices.empty()) 
            {
                size_t rowIdx = 0;
                for (size_t i : promisingIndices) 
                {
                    const auto& pt = historyPoints[i];
                    for (size_t j = 0; j < numNeighbors; ++j) 
                    {
                        for (size_t d = 0; d < N; ++d) 
                        {
                            double shift = ((j >> d) & 1) ? delta[d] : -delta[d];
                            double neighborVal = pt[d] + shift;
                            neighborVal = std::max(static_cast<double>(lowerBound[d]), std::min(static_cast<double>(upperBound[d]), neighborVal));
                            
                            neighborSamplesBuffer.at<float>(rowIdx, d) = static_cast<float>(neighborVal);
                        }
                        rowIdx++;
                    }
                }

                cv::Mat activeSamples = neighborSamplesBuffer.rowRange(0, static_cast<int>(rowIdx));
                randomForest->predict(activeSamples, predictions);

                for (size_t idx = 0; idx < promisingIndices.size(); ++idx) 
                {
                    size_t ptIndex = promisingIndices[idx];
                    double ptVal = historyValues[ptIndex];
                    bool isLocalExtremum = true;
                    
                    size_t baseRow = idx * numNeighbors;
                    for (size_t j = 0; j < numNeighbors; ++j) 
                    {
                        double predZ = predictions.at<float>(baseRow + j, 0);
                        if (isMinimize ? (ptVal >= predZ) : (ptVal <= predZ)) 
                        { 
                            isLocalExtremum = false; 
                            break; 
                        }
                    }

                    if (isLocalExtremum) 
                    {
                        LocalMethod<N> lm(func, historyPoints[ptIndex], lowerBound, upperBound, isMinimize);
                        lm.SetMaxTrial(40); 
                        auto localRes = lm.StartOptimization();
                        
                        double internalLocalZ = isMinimize ? static_cast<double>(localRes.second) : -static_cast<double>(localRes.second);

                        if (!foundAnyCandidate || 
                           (isMinimize ? (internalLocalZ < overallBestLocalZ) : (internalLocalZ > overallBestLocalZ))) 
                        {
                            overallBestLocalZ = internalLocalZ;
                            overallBestLocalPoint = localRes.first;
                            foundAnyCandidate = true;
                        }
                    }
                }
            }

            if (foundAnyCandidate) 
            {
                int xEvalLocalIndex = 0;
                
                for (size_t d = 0; d < N; d++) 
                {
                    bestLocalDouble[d] = overallBestLocalPoint[d];
                }
                
                double trueZLocalInternal = sign * func->ComputePoint(bestLocalDouble, xEvalLocalIndex);

                if (isMinimize ? (trueZLocalInternal < resZInternal) : (trueZLocalInternal > resZInternal)) 
                {
                    if (xEvalLocalIndex == evalIndex) 
                    {
                        resZInternal = trueZLocalInternal;
                        bestPoint = overallBestLocalPoint;
                    }
                }

                if (xEvalLocalIndex < Z.size() && trueZLocalInternal < Z[xEvalLocalIndex])
                {
                    Z[xEvalLocalIndex] = trueZLocalInternal;
                }

                historyPoints.push_back(overallBestLocalPoint);
                historyValues.push_back(trueZLocalInternal);
            }
        }

    }   while ((std::abs(resZInternal - resPoint) > eps) && (iteration < maxIteration));

    double finalZ = isMinimize ? resZInternal : -resZInternal;
    std::vector<double> result;
    result.push_back(finalZ);
    for (size_t i = 0; i < N; ++i) 
    {
        result.push_back(bestPoint[i]);
    }
    result.push_back(static_cast<double>(iteration));
    
    return result;
}