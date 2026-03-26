#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>

#include "Algorithm.h"
#include "HillProblem.hpp"
#include "ShekelProblem.hpp"
#include "grishagin_function.hpp"
#include "GKLSProblem.hpp"
#include "GKLSConstrainedProblem.hpp"

// =========================================================================
// Универсальный раннер для статических семейств (Hill, Shekel, Grishagin)
// =========================================================================
template<typename ProblemClass, size_t N>
void RunFamilyTest(const std::string& name, int startIndex, int count, double r, int tightness, double eps, bool isMinimize) 
{
    size_t success = 0;
    size_t fails = 0;
    size_t totalIters = 0;
    
    std::cout << std::left << std::setw(15) << name 
              << " | " << std::setw(3) << (isMinimize ? "MIN" : "MAX")
              << " | r=" << std::fixed << std::setprecision(1) << std::setw(3) << r << " | ";

    for (int i = startIndex; i < startIndex + count; ++i) 
    {
        ProblemClass prob(i);
        double expZ = 0.0;
        std::vector<double> expX;
        
        try 
        {
            if (isMinimize) {
                if (!prob.GetStatus(ofpMinimum)) throw std::runtime_error("No min");
                expZ = prob.GetOptimumValue();
                expX = prob.GetOptimumPoint();
            } else {
                if (!prob.GetStatus(ofpMaximum)) throw std::runtime_error("No max");
                expZ = prob.GetMaxValue();
                expX = prob.GetMaxPoint();
            }
        } 
        catch (...) 
        {
            std::cout << "Экстремум не задан (тест пропущен)\n";
            return; 
        }

        std::vector<double> lb, ub;
        prob.GetBounds(lb, ub);
        TPoint<double, N> lower(lb);
        TPoint<double, N> upper(ub);
        
        TAlgorithm<double, N> alg(lower, upper, eps, r, &prob, tightness);
        
        size_t maxIters = (N == 1) ? 5000 : 40000;
        std::vector<double> res = alg.Solve(maxIters, isMinimize);
        
        if (res.empty() || res.size() < N + 2) {
            fails++;
            continue;
        }
        
        double gotZ = res[0];
        bool pointOk = false;
        
        if (expX.size() >= N) {
            double maxDist = 0.0;
            for (size_t d = 0; d < N; ++d) {
                maxDist = std::max(maxDist, std::abs(res[d+1] - expX[d]));
            }
            pointOk = (maxDist <= 0.05);
        }
        
        bool valOk = (std::abs(gotZ - expZ) <= 0.05);
        
        if (pointOk || valOk) {
            success++;
            totalIters += static_cast<size_t>(res.back());
        } else {
            fails++;
        }
    }
    
    double rate = (double)success / count * 100.0;
    double avgIter = success > 0 ? (double)totalIters / success : 0;
    
    std::cout << "Успех: " << std::setw(6) << std::setprecision(2) << rate << "% | "
              << "Ср. итераций: " << std::setw(5) << std::setprecision(0) << avgIter << "\n";
}

// =========================================================================
// Раннер для безусловных задач GKLS (у них своя сигнатура конструктора)
// =========================================================================
template<size_t N>
void RunGKLSTest(const std::string& name, GKLSClass gklsClass, int count, double r, int tightness, double gkls_eps) 
{
    size_t success = 0, fails = 0, totalIters = 0;
    
    std::cout << std::left << std::setw(15) << name 
              << " | MIN | r=" << std::fixed << std::setprecision(1) << std::setw(3) << r << " | ";

    // Пересчет допуска для многомерной развертки Пеано
    double alg_eps = std::pow(gkls_eps / 4.0, static_cast<double>(N)); 

    for (int i = 1; i <= count; ++i) 
    {
        TGKLSProblem prob(i, static_cast<int>(N), gklsClass, TD);
        
        std::vector<double> lb, ub;
        prob.GetBounds(lb, ub);
        TPoint<double, N> lower(lb);
        TPoint<double, N> upper(ub);
        
        TAlgorithm<double, N> alg(lower, upper, alg_eps, r, &prob, tightness);
        
        size_t maxIters = (N >= 4) ? 2000000 : 500000;
        std::vector<double> res = alg.Solve(maxIters, true);
        
        if (res.empty() || res.size() < N + 2) {
            fails++;
            continue;
        }
        
        std::vector<double> expX = prob.GetOptimumPoint();
        double gotZ = res[0];
        double expZ = prob.GetOptimumValue();
        
        bool pointOk = false;
        if (expX.size() >= N) {
            double maxDist = 0.0;
            for (size_t d = 0; d < N; ++d) {
                maxDist = std::max(maxDist, std::abs(res[d+1] - expX[d]));
            }
            pointOk = (maxDist <= gkls_eps);
        }
        
        bool valOk = (std::abs(gotZ - expZ) <= gkls_eps);
        
        if (pointOk || valOk) {
            success++;
            totalIters += static_cast<size_t>(res.back());
        } else {
            fails++;
        }
    }
    
    double rate = (double)success / count * 100.0;
    double avgIter = success > 0 ? (double)totalIters / success : 0;
    
    std::cout << "Успех: " << std::setw(6) << std::setprecision(2) << rate << "% | "
              << "Ср. итераций: " << std::setw(5) << std::setprecision(0) << avgIter << "\n";
}

inline void RunConsoleBenchmarks() 
{
    std::cout << "===================================================================\n";
    std::cout << "   СВОДНОЕ ТЕСТИРОВАНИЕ: HILL, SHEKEL, GRISHAGIN, GKLS \n";
    std::cout << "===================================================================\n";
    
    RunGKLSTest<3>("GKLS 3D Simple", Simple, 100, 4.5, 7,  0.01);
    RunGKLSTest<3>("GKLS 3D Hard",   Hard,   100, 5.6, 7,  0.01);

    std::cout << "===================================================================\n";
}