#ifndef __LOCALMETHOD_H__
#define __LOCALMETHOD_H__

#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>

template<size_t N>
class LocalMethod
{
protected:
    IGeneralOptProblem* mFunc;
    int mConstraintsNumber;
    int mTrialsCounter;
    int mMaxTrial;

    TPoint<N> mBestPoint;
    double mBestValue;
    
    double mEps;
    double mStep;
    double mStepMultiplier;

    double mStartPoint[N];
    double mCurrentPoint[N];
    double mCurrentResearchDirection[N];
    double mPreviousResearchDirection[N];
    
    TPoint<N> mLowerBound;
    TPoint<N> mUpperBound;
    bool mIsMinimize;

    double EvaluateObjectiveFunction(const double* x)
    {
        if (mTrialsCounter >= mMaxTrial)
        {
            return HUGE_VAL;
        }

        for (size_t i = 0; i < N; i++) 
        {
            if (x[i] < static_cast<double>(mLowerBound[i]) || x[i] > static_cast<double>(mUpperBound[i]))
            {
                return HUGE_VAL;
            }
        }

        int va = 0;
        std::vector<double> point(x, x + N); 
        double value = mFunc->ComputePoint(point, va);
        mTrialsCounter++;
        
        return mIsMinimize ? value : -value; 
    }

    double MakeResearch(double* startPoint)
    {
        double currentFValue = EvaluateObjectiveFunction(startPoint);
        for (size_t i = 0; i < N; i++) 
        {
            startPoint[i] += mStep;
            double rightFValue = EvaluateObjectiveFunction(startPoint);
            
            if (rightFValue < currentFValue) 
            {
                currentFValue = rightFValue;
            }
            else 
            {
                startPoint[i] -= 2 * mStep;
                double leftFValue = EvaluateObjectiveFunction(startPoint);
                
                if (leftFValue < currentFValue) 
                {
                    currentFValue = leftFValue;
                }
                else 
                {
                    startPoint[i] += mStep; 
                }
            }
        }
        return currentFValue;
    }

    void DoStep()
    {
        for (size_t i = 0; i < N; i++) 
        {
            mCurrentPoint[i] = mCurrentResearchDirection[i] + (mCurrentResearchDirection[i] - mPreviousResearchDirection[i]);
        }
    }

public:
    LocalMethod(IGeneralOptProblem* func, TPoint<N> startPoint, TPoint<N> lb, TPoint<N> ub, bool isMinimize)
        : mFunc(func), mBestPoint(startPoint), mLowerBound(lb), mUpperBound(ub), mIsMinimize(isMinimize)
    {
        mConstraintsNumber = 0;
        mTrialsCounter = 0;
        mMaxTrial = 1000;
        mEps = 0.001;
        mStep = 0.1;
        mStepMultiplier = 2.0;

        for (size_t i = 0; i < N; ++i) 
        {
            mStartPoint[i] = static_cast<double>(startPoint[i]);
        }

        mBestValue = EvaluateObjectiveFunction(mStartPoint);
    }

    void SetMaxTrial(int maxTrial) 
    {
        mMaxTrial = maxTrial;
    }

    std::pair<TPoint<N>, double> StartOptimization()
    {
        int i = 0;
        int k = 0;
        bool needRestart = true;
        double currentFValue = mBestValue;
        mTrialsCounter = 0;

        while (mTrialsCounter < mMaxTrial) 
        {
            i++;
            if (needRestart) 
            {
                k = 0;
                std::memcpy(mCurrentPoint, mStartPoint, sizeof(double) * N);
                std::memcpy(mCurrentResearchDirection, mStartPoint, sizeof(double) * N);
                currentFValue = EvaluateObjectiveFunction(mCurrentPoint);
                needRestart = false;
            }

            std::swap(mPreviousResearchDirection, mCurrentResearchDirection);
            std::memcpy(mCurrentResearchDirection, mCurrentPoint, sizeof(double) * N);
            double nextFValue = MakeResearch(mCurrentResearchDirection);

            if (currentFValue > nextFValue) 
            {
                DoStep();
                k++;
                currentFValue = nextFValue;
            }
            else if (mStep > mEps) 
            {
                if (k != 0)
                {
                    std::memcpy(mStartPoint, mPreviousResearchDirection, sizeof(double) * N);
                }
                else
                {
                    mStep /= mStepMultiplier;
                }
                
                needRestart = true;
            }
            else
            {
                break;
            }
        }

        if (currentFValue < mBestValue) 
        {
            for (size_t j = 0; j < N; ++j) 
            {
                mBestPoint[j] = static_cast<double>(mPreviousResearchDirection[j]);
            }
            mBestValue = currentFValue;
        }

        return {mBestPoint, mIsMinimize ? mBestValue : -mBestValue};
    }
};

#endif