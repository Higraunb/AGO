#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "Algorithm.h"
#include "HillProblem.hpp"
#include "ShekelProblem.hpp"

using namespace std;

template <typename ProblemClass>
void run_tests(double testarr[][2], 
               double lower, double upper,
               double eps, size_t max_tests, 
               const string& testName, const double r)
{
    size_t countSuccess = 0;
    size_t countFail = 0;
    vector<vector<double>> failTests;
    
    TPoint<double, 1> lowerBoundVec; lowerBoundVec[0] = lower;
    TPoint<double, 1> upperBoundVec; upperBoundVec[0] = upper;

    cout << "=======================================\n";
    cout << "STARTING TESTS: " << testName << "\n";
    cout << "Range: [" << lower << ", " << upper << "]\n";
    cout << "=======================================\n";

    for (size_t i = 0; i < max_tests; i++)
    {
        ProblemClass func(i);

        TAlgorithm<double, 1> alg(lowerBoundVec, upperBoundVec, eps, r, &func);
        vector<double> res = alg.AGPStronginaMin();

        double expectedX = testarr[i][0];
        double expectedZ = testarr[i][1];

        bool xOk = abs(expectedX - res[0]) < 0.05;
        bool zOk = abs(expectedZ - res[1]) < 0.05;
        if (zOk && xOk)
        {
            countSuccess++;
            // cout << "Test " << i << ": SUCCESS (iter: " << res[2] << ")\n";
        }
        else
        {
            countFail++;
            failTests.push_back({ (double)i, expectedX, res[0], expectedZ, res[1], res[2] });
            cout << "Test " << i << ": FAIL\n";
        }
    }

    cout << "\nFailed tests summary (" << testName << "):\n";
    for (const auto& f : failTests)
    {
        cout << "Test " << (int)f[0] 
             << " | Exp X: " << f[1] << ", Got X: " << f[2] 
             << " | Exp Z: " << f[3] << ", Got Z: " << f[4] 
             << " | Iters: " << f[5] << "\n";
    }

    cout << "---------------------------------------\n";
    cout << "Total: " << max_tests << " | Success: " << countSuccess << " | Fail: " << countFail << "\n\n";
}

int main() 
{
    double eps = 0.0001;

    // run_tests<THillProblem>(minHill, 0.0, 1.0, eps, 1000, "HILL Problem", 4); // r = 4.0

    run_tests<TShekelProblem>(minShekel, 0.0, 10.0, eps, 1000, "SHEKEL Problem", 7); // r = 7.0

    return 0;
}