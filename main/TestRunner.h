#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <locale> 
#include "Algorithm.h"

#include "HillProblem.hpp"
#include "ShekelProblem.hpp"

using namespace std;

struct CommaNumPunct : std::numpunct<char> {
    char do_decimal_point() const override { return ','; }
};

struct TestStats {
    size_t total;
    size_t success;
    size_t fail;
    double successRate;
    double avgIterations;
};

struct FailedTest {
    size_t id;
    double expZ;
    double gotZ;
    double expX;
    double gotX;
};

template <typename ProblemClass>
TestStats run_single_r_test(double testarr[][2], double lower, double upper, double eps, size_t max_tests, double r, int tightness_)
{
    size_t countSuccess = 0;
    size_t countFail = 0;
    size_t totalIterations = 0;
    
    vector<FailedTest> failures;
    
    TPoint<double, 1> lowerBoundVec; lowerBoundVec[0] = lower;
    TPoint<double, 1> upperBoundVec; upperBoundVec[0] = upper;

    for (size_t i = 1; i < max_tests; i++)
    {
        ProblemClass func(i);
        TAlgorithm<double, 1> alg(lowerBoundVec, upperBoundVec, eps, r, &func, tightness_);
        vector<double> res = alg.Solve(1000, true); 

        double expectedX = testarr[i][1];
        double expectedZ = testarr[i][0];

        bool xOk = abs(expectedX - res[1]) < 0.1;
        bool zOk = abs(expectedZ - res[0]) < 0.1;
        
        if (zOk && xOk)
        {
            countSuccess++;
            totalIterations += (size_t)res[2];
        }
        else
        {
            countFail++;
            failures.push_back({i, expectedZ, res[0], expectedX, res[1]});
        }
    }


    if (countFail == 1 || countFail == 2)
    {
        cout << "  [!] Внимание: при r = " << fixed << setprecision(1) << r << " упало тестов: " << countFail << "\n";
        for (const auto& f : failures)
        {
            cout << "      Тест " << setw(3) << f.id 
                 << " | Ожидали X: " << setw(8) << f.expX << " | Получили X: " << setw(8) << f.gotX
                 << " | Ожидали Z: " << setw(8) << f.expZ << " | Получили Z: " << setw(8) << f.gotZ << "\n";
        }
    }

    double successRate = (double)countSuccess / (max_tests - 1) * 100.0;
    double avgIter = countSuccess > 0 ? (double)totalIterations / countSuccess : 0.0;

    return {max_tests - 1, countSuccess, countFail, successRate, avgIter};
}

template <typename ProblemClass>
void run_experiment(double testarr[][2], double lower, double upper, double eps, size_t max_tests, const string& testName, const string& filename)
{
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Ошибка: не удалось открыть файл " << filename << " для записи!\n";
        return;
    }


    file.imbue(std::locale(file.getloc(), new CommaNumPunct()));

    file << "r;Success_Rate_%;Success_Count;Fail_Count;Avg_Iterations\n";

    cout << "============================================================\n";
    cout << "  STARTING EXPERIMENT: " << testName << "\n";
    cout << "  Saving results to: " << filename << "\n";
    cout << "============================================================\n";

    for (double r = 1.5; r <= 8.01; r += 0.2) 
    {
        int tightness = 3;
        TestStats stats = run_single_r_test<ProblemClass>(testarr, lower, upper, eps, max_tests, r, tightness);
        
        file << r << ";" 
             << stats.successRate << ";" 
             << stats.success << ";" 
             << stats.fail << ";" 
             << stats.avgIterations << "\n";

        cout << "Tested r = " << fixed << setprecision(1) << setw(3) << r 
             << " | Rate: " << setw(6) << setprecision(2) << stats.successRate << "%"
             << " | Avg Iters: " << setw(6) << setprecision(1) << stats.avgIterations << "\n";
    }
    cout << "============================================================\n";
    cout << "Эксперимент завершен. Данные сохранены в " << filename << "\n\n";
    
    file.close();
}

inline void runAllTests() 
{
    double eps = 0.01;
    
    run_experiment<TShekelProblem>(minShekel, 0.0, 10.0, eps, 1000, "SHEKEL Problem", "Shekel_Results.csv");
}