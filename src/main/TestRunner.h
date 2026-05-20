#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <chrono>

#include "Algorithm.h"
#include "HillProblem.hpp"
#include "ShekelProblem.hpp"
#include "grishagin_function.hpp"
#include "GKLSProblem.hpp"
#include "GKLSConstrainedProblem.hpp"
#include "../../../Logger/Logger.h" 

struct TestStats 
{
    size_t total;
    size_t success;
    size_t fail;
    double successRate;
    double avgIterations;
    double totalTimeSec;
    double avgTimeSec;
};

template<typename ProblemClass, size_t N>
void RunFamilyTest(const std::string& name, int startIndex, int count, double r, int tightness, double eps, bool isMinimize) 
{
    size_t success = 0;
    size_t fails = 0;
    size_t totalIters = 0;

    for (int i = startIndex; i < startIndex + count; ++i) 
    {
        ProblemClass prob(i);
        double expZ = 0.0;
        std::vector<double> expX;
        
        try 
        {
            if (isMinimize) 
            {
                if (!prob.GetStatus(ofpMinimum)) throw std::runtime_error("No min");
                expZ = prob.GetOptimumValue();
                expX = prob.GetOptimumPoint();
            } 
            else 
            {
                if (!prob.GetStatus(ofpMaximum)) throw std::runtime_error("No max");
                expZ = prob.GetMaxValue();
                expX = prob.GetMaxPoint();
            }
        } 
        catch (...) 
        {
            continue; 
        }

        std::vector<double> lb, ub;
        prob.GetBounds(lb, ub);
        TPoint<N> lower(lb);
        TPoint<N> upper(ub);
        
        Algorithm<N> alg(lower, upper, eps, r, &prob, tightness);
        
        size_t maxIters = (N == 1) ? 5000 : 40000;
        int p = 150;
        std::vector<double> res = alg.Solve(maxIters, p, isMinimize);
          
        if (res.empty() || res.size() < N + 2) 
        {
            fails++;
            LOG_ERROR("  [ERROR] {} task #{}: Algorithm returned empty result!", name, i);
            continue;
        }
        
        double gotZ = res[0];
        bool pointOk = false;
        double maxDist = 0.0;
        
        if (expX.size() >= N) 
        {
            for (size_t d = 0; d < N; ++d) 
            {
                maxDist = std::max(maxDist, std::abs(res[d+1] - expX[d]));
            }
            pointOk = (maxDist <= 0.05);
        }
        
        bool valOk = (std::abs(gotZ - expZ) <= 0.05);
        
        if (pointOk || valOk) 
        {
            success++;
            totalIters += static_cast<size_t>(res.back());
        } 
        else 
        {
            fails++;
            LOG_INFO("  [MISS] {} task #{}: Expected Z={:.6f}, Found Z={:.6f} | Diff X={:.6f} (tolerance 0.05)", 
                     name, i, expZ, gotZ, maxDist);
        }
    }
    
    double rate = (double)success / count * 100.0;
    double avgIter = success > 0 ? (double)totalIters / success : 0;
    
    std::string result_msg = fmt::format("{:<15} | {:<3} | r={:<3.1f} | Success: {:>6.2f}% | Avg. iterations: {:>5.0f}", 
                                         name, (isMinimize ? "MIN" : "MAX"), r, rate, avgIter);
    LOG_INFO(result_msg);
}

template<size_t N>
TestStats RunGKLSTest(const std::string& name, GKLSClass gklsClass, int count, double r, int tightness, double gkls_eps) 
{
    size_t success = 0, fails = 0, totalIters = 0;

    double alg_eps = gkls_eps;
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 1; i <= count; ++i) 
    {
        TGKLSProblem prob(i, static_cast<int>(N), gklsClass, TD);
        
        std::vector<double> lb, ub;
        prob.GetBounds(lb, ub);
        TPoint<N> lower(lb);
        TPoint<N> upper(ub);
        
        Algorithm<N> alg(lower, upper, alg_eps, r, &prob, tightness);
        
        size_t maxIters = (N >= 4) ? 2000000 : 500000;
        int p = 150;
        std::vector<double> res = alg.Solve(maxIters, p, true);
        
        if (res.empty() || res.size() < N + 2) 
        {
            fails++;
            LOG_ERROR("  [ERROR] {} task #{}: Algorithm returned empty or incomplete result!", name, i);
            continue;
        }
        
        std::vector<double> expX = prob.GetOptimumPoint();
        double gotZ = res[0];
        double expZ = prob.GetOptimumValue();
        
        bool pointOk = false;
        double maxDist = 0.0;
        if (expX.size() >= N) 
        {
            for (size_t d = 0; d < N; ++d) 
            {
                maxDist = std::max(maxDist, std::abs(res[d+1] - expX[d]));
            }
            pointOk = (maxDist <= gkls_eps);
        }
        
        bool valOk = (std::abs(gotZ - expZ) <= gkls_eps);
        
        if (pointOk || valOk) 
        {
            success++;
            totalIters += static_cast<size_t>(res.back());
        } 
        else 
        {
            fails++;
            LOG_INFO("  [MISS] {} task #{}: Expected Z={:.6f}, Found Z={:.6f} | X={:.6f} Y={:.6f} (tolerance {:.6f}) iteration {:.1f}", 
                     name, i, expZ, gotZ, res[1], res[2], gkls_eps, res[3]);
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;
    double totalTime = diff.count();
    double avgTime = totalTime / count;

    double rate = (double)success / count * 100.0;
    double avgIter = success > 0 ? (double)totalIters / success : 0;
    
    std::string result_msg = fmt::format("{:<15} | MIN | r={:<3.1f} | Success: {:>6.2f}% | Avg. iterations: {:>5.0f} | Time: {:.2f}s", 
                                         name, r, rate, avgIter, totalTime);
    LOG_INFO(result_msg);

    return { (size_t)count, success, fails, rate, avgIter, totalTime, avgTime };
}

inline void RunConsoleBenchmarks() 
{
    double eps = 0.001; 
    int tightness = 10; 

    LOG_INFO("===================================================================");
    LOG_INFO("   GKLS SUMMARY TESTING");
    LOG_INFO("===================================================================");
    
    std::vector<int>    dims    = { 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 8 };
    std::vector<double> epss    = { 0.002, 0.002, 0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.0464, 0.0517, 0.0562 };
    std::vector<std::string> shs     = { "simple", "hard", "simple", "hard", "simple", "hard", "simple", "hard", "simple", "simple", "simple" };
    std::vector<double> rs      = { 5, 6.7, 4.5, 5.6, 4.5, 4.6, 4.5, 5.6, 5.6, 5.6, 5.6 };
    
    for (int i = 0; i <= 7; i++) 
    {
        int N = dims[i];
        double gkls_eps = epss[i];
        std::string sh0 = shs[i];
        double r = rs[i];
        
        GKLSClass type = (sh0 == "simple") ? Simple : Hard;
        std::string name = fmt::format("GKLS {}D {}", N, sh0 == "simple" ? "Simple" : "Hard  ");

        switch (N) 
        {
            case 2: RunGKLSTest<2>(name, type, 100, r, tightness, gkls_eps); break;
            case 3: RunGKLSTest<3>(name, type, 100, r, tightness, gkls_eps); break;
            case 4: RunGKLSTest<4>(name, type, 100, r, tightness, gkls_eps); break;
            case 5: RunGKLSTest<5>(name, type, 100, r, tightness, gkls_eps); break;
            default: LOG_ERROR("Dimension {} is not supported in switch", N); continue;
        }
    }
    
    LOG_INFO("===================================================================");
    LOG_INFO("   Experiments completed.");
    LOG_INFO("===================================================================");
}