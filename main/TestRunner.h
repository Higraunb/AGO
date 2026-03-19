#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <locale>
#include <cmath>
#include "Algorithm.h"
#include "grishagin_function.hpp"
#include <filesystem>
#include "GKLSProblem.hpp"
#include "GKLSConstrainedProblem.hpp"

using namespace std;

// Локаль для Excel (чтобы выводились запятые)
struct CommaNumPunct : std::numpunct<char> {
    char do_decimal_point() const override { return ','; }
};

struct TestStats {
    size_t total;
    size_t success;
    size_t fail;
    double successRate;
    double avgIterations;
    double totalTimeSec;
    double avgTimeSec;
};

// Обновленная структура для хранения данных упавшего теста (теперь есть X и Y)
struct FailedTest {
    size_t id;
    double expZ;
    double gotZ;
    double expX;
    double gotX;
    double expY;
    double gotY;
};

template <size_t N>
TestStats run_constrained_batch(double gkls_eps, double r, GKLSClass type,
                                size_t max_tests, int tightness, int constrType)
{
    size_t countSuccess    = 0;
    size_t countFail       = 0;
    size_t totalIterations = 0;

    TPoint<double, N> lowerBoundVec;
    TPoint<double, N> upperBoundVec;
    for (size_t i = 0; i < N; ++i) {
        lowerBoundVec[i] = -1.0;
        upperBoundVec[i] =  1.0;
    }

    double alg_eps = std::pow(gkls_eps / 4.0, (int)N);

    auto start_time = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for reduction(+:countSuccess, countFail, totalIterations)
    for (size_t i = 1; i <= max_tests; i++)
    {
        TGKLSConstrainedProblem func(
            static_cast<EConstrainedProblemType>(constrType),
            0.5, 0, (int)i, (int)N, type, TD);

        TAlgorithm<double, N> alg(lowerBoundVec, upperBoundVec, alg_eps, r, &func, tightness);

        size_t max_iters = (N >= 4) ? 2000000 : 500000;
        vector<double> res = alg.Solve(max_iters, true);

        vector<double> foundPoint(res.begin() + 1, res.begin() + 1 + N);

        // --- Критерий 1: точка лежит в допустимой области ---
        int finalIndex = 0;
        func.ComputePoint(foundPoint, finalIndex);
        bool isFeasible = (finalIndex == func.GetConstraintsNumber());

        // --- Критерий 2: близость к известному глобальному оптимуму ---
        // GetOptimumPoint() бросает исключение для cptNormal если точка нарушает
        // сырые ограничения — оборачиваем в try/catch
        bool isGlobalOptimum = false;
        try
        {
            vector<double> opt_point = func.GetOptimumPoint();
            double max_dist = 0.0;
            for (size_t d = 0; d < N; ++d)
                max_dist = std::max(max_dist,
                                    std::abs(foundPoint[d] - opt_point[d]));
            isGlobalOptimum = (max_dist <= gkls_eps);
        }
        catch (const std::string&)
        {
            // GetOptimumPoint недоступен для данного типа задачи —
            // считаем успех только по feasibility + значению целевой функции
            if (isFeasible)
            {
                // Сравниваем найденное значение с оптимальным
                double foundZ  = res[0];
                double optZ    = func.GetOptimumValue();
                isGlobalOptimum = (std::abs(foundZ - optZ) <= gkls_eps);
            }
        }

        if (isFeasible && isGlobalOptimum) {
            countSuccess++;
            totalIterations += (size_t)res[N + 1];
        } else {
            countFail++;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;
    double totalTime = diff.count();
    double avgTime   = totalTime / max_tests;

    double successRate = (double)countSuccess / max_tests * 100.0;
    double avgIter     = countSuccess > 0
                         ? (double)totalIterations / countSuccess : 0.0;

    return {max_tests, countSuccess, countFail,
            successRate, avgIter, totalTime, avgTime};
}

template <typename ProblemClass>
TestStats run_single_r_test_2d(double lower, double upper, double eps, size_t max_tests, double r, int tightness)
{
    size_t countSuccess = 0;
    size_t countFail = 0;
    size_t totalIterations = 0;
    
    vector<FailedTest> failures;
    
    // Для задачи Гришагина размерность N = 2
    TPoint<double, 2> lowerBoundVec(lower, lower);
    TPoint<double, 2> upperBoundVec(upper, upper);

    // Задачи Гришагина нумеруются строго от 1 до 100 включительно
    for (size_t i = 1; i <= max_tests; i++)
    {
        ProblemClass func(i);
        
        // Передаем N = 2
        TAlgorithm<double, 2> alg(lowerBoundVec, upperBoundVec, eps, r, &func, tightness);
        
        // Запускаем на 10000 итераций (многомерные задачи требуют больше шагов)
        vector<double> res = alg.Solve(40000, true); 

        // Достаем ожидаемые координаты X и Y из глобального массива Гришагина
        // Индекс массива начинается с 0, поэтому (i - 1)
        double expectedX = rand_minimums[2 * (i - 1)];
        double expectedY = rand_minimums[2 * (i - 1) + 1];
        
        // Вычисляем ожидаемое значение функции в точке глобального минимума
        std::vector<double> expPoint = { expectedX, expectedY };
        double expectedZ = func.ComputeFunction(expPoint);

        // Разбираем ответ от алгоритма. При N=2: [0]-X, [1]-Y, [2]-Z, [3]-Iters
        double gotZ = res[0];
        double gotX = res[1];
        double gotY = res[2];
        size_t iters = (size_t)res[3];

        // Допуск для многомерной оптимизации с разверткой (эвольвентой) обычно делают чуть шире
        bool xOk = abs(expectedX - gotX) < 0.05;
        bool yOk = abs(expectedY - gotY) < 0.05;
        
        if (xOk && yOk)
        {
            countSuccess++;
            totalIterations += iters;
        }
        else
        {
            countFail++;
            failures.push_back({i, expectedZ, gotZ, expectedX, gotX, expectedY, gotY});
        }
    }

    // Выводим информацию, если упало малое количество тестов (например, до 3)
    /*if (countFail > 0 && countFail <= 3)
    {
        cout << "  [!] Внимание: при r = " << fixed << setprecision(1) << r << " упало тестов: " << countFail << "\n";
        for (const auto& f : failures)
        {
            cout << "      Тест " << setw(3) << f.id 
                 << " | Ожидали: (" << setw(5) << f.expX << ", " << setw(5) << f.expY << ")"
                 << " | Получили: (" << setw(5) << f.gotX << ", " << setw(5) << f.gotY << ")"
                 << " | Z_Ожид: " << setw(8) << f.expZ << " | Z_Получ: " << setw(8) << f.gotZ << "\n";
        }
    }*/

    double successRate = (double)countSuccess / max_tests * 100.0;
    double avgIter = countSuccess > 0 ? (double)totalIterations / countSuccess : 0.0;

    return {max_tests, countSuccess, countFail, successRate, avgIter};
}


inline void runConstrainedTests()
{
    vector<int>    dims = { 2, 3, 4 };
    vector<double> epss = { 0.01, 0.01, 0.03 };
    vector<string> shs  = { "simple", "hard", "simple" };
    vector<double> rs   = { 4.5, 5.6, 4.5 };

    // Исправлено: 1=cptInFeasibleDomain, 3=cptOnFeasibleBorder
    vector<int>    constrTypes = { 1, 3 };
    vector<string> constrNames = { "InFeasibleDomain", "OnFeasibleBorder" };

    std::filesystem::create_directories("results");
    string filename = "results/Constrained_GKLS_Results.csv";

    ofstream file(filename);
    if (file.is_open()) {
        file.imbue(std::locale(file.getloc(), new CommaNumPunct()));
        file << "Type;Dim;Class;eps;r;Success_Rate_%;Success_Count;"
                "Fail_Count;Avg_Iterations;Total_Time_s;Avg_Time_s\n";
    }

    cout << "============================================================\n";
    cout << "  STARTING CONSTRAINED EXPERIMENTS\n";
    cout << "============================================================\n";

    for (size_t c = 0; c < constrTypes.size(); c++)
    {
        cout << "--- Constraint Type: " << constrNames[c] << " ---\n";

        for (size_t i = 0; i < dims.size(); i++)
        {
            int    N        = dims[i];
            double gkls_eps = epss[i];
            string sh0      = shs[i];
            double r        = rs[i];

            GKLSClass type  = (sh0 == "simple") ? Simple : Hard;

            int tightness = std::ceil(-std::log2(gkls_eps / 4.0));
            if (tightness * N > 53) tightness = 53 / N;

            cout << "C_GKLS | Dim=" << N << " | " << setw(6) << sh0
                 << " | r=" << fixed << setprecision(1) << r
                 << " | tight=" << tightness << " ... " << flush;

            TestStats stats{};

            if      (N == 2) stats = run_constrained_batch<2>(gkls_eps, r, type, 100, tightness, constrTypes[c]);
            else if (N == 3) stats = run_constrained_batch<3>(gkls_eps, r, type, 100, tightness, constrTypes[c]);
            else if (N == 4) stats = run_constrained_batch<4>(gkls_eps, r, type, 100, tightness, constrTypes[c]);

            if (file.is_open()) {
                file << constrNames[c] << ";" << N << ";" << sh0 << ";"
                     << gkls_eps << ";" << r << ";"
                     << stats.successRate    << ";" << stats.success << ";"
                     << stats.fail           << ";" << stats.avgIterations << ";"
                     << stats.totalTimeSec   << ";" << stats.avgTimeSec << "\n";
            }

            cout << (stats.successRate == 100.0 ? "PERFECT!" : "Done.")
                 << " Rate: " << setw(6) << setprecision(2) << stats.successRate
                 << "% | Iters: " << fixed << setprecision(0) << stats.avgIterations
                 << " | Time: " << setprecision(2) << stats.totalTimeSec << "s\n";
        }
    }

    cout << "============================================================\n";
    cout << "Эксперименты завершены. Результаты в " << filename << "\n\n";
    if (file.is_open()) file.close();
}

template <typename ProblemClass>
void run_experiment_2d(double lower, double upper, double eps, size_t max_tests, const string& testName, const string& filename)
{
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Ошибка: не удалось открыть файл " << filename << " для записи!\n";
        return;
    }

    file.imbue(std::locale(file.getloc(), new CommaNumPunct()));
    file << "r;Success_Rate_%;Success_Count;Fail_Count;Avg_Iterations\n";

    cout << "============================================================\n";
    cout << "  STARTING EXPERIMENT: " << testName << " (2D Dimension)\n";
    cout << "  Saving results to: " << filename << "\n";
    cout << "============================================================\n";

    // Плотность сетки (эвольвенты) Пеано для 2D задач. Чем больше - тем точнее, но медленнее.
    int tightness = 7;

    for (double r = 2.0; r <= 10.0; r += 0.2) 
    {
        TestStats stats = run_single_r_test_2d<ProblemClass>(lower, upper, eps, max_tests, r, tightness);
        
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
template <size_t N>
TestStats run_gkls_batch(double gkls_eps, double r, GKLSClass type, size_t max_tests, int tightness)
{
    size_t countSuccess = 0;
    size_t countFail = 0;
    size_t totalIterations = 0;

    TPoint<double, N> lowerBoundVec;
    TPoint<double, N> upperBoundVec;
    for (size_t i = 0; i < N; ++i) {
        lowerBoundVec[i] = -1.0;
        upperBoundVec[i] = 1.0;
    }

    double alg_eps = std::pow(gkls_eps / 4.0, N); 

    // Запускаем таймер высокого разрешения перед началом пачки тестов
    auto start_time = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for reduction(+:countSuccess, countFail, totalIterations)
    for (size_t i = 1; i <= max_tests; i++)
    {
        TGKLSProblem func(i, N, type, TD);
        TAlgorithm<double, N> alg(lowerBoundVec, upperBoundVec, alg_eps, r, &func, tightness);
        
        size_t max_iters = (N >= 4) ? 2000000 : 500000; 
        vector<double> res = alg.Solve(max_iters, true); 

        vector<double> opt_point = func.GetOptimumPoint();
        
        bool isOk = true;
        double max_dist = 0.0;
        for (size_t d = 0; d < N; ++d) {
            double dist = std::abs(res[d + 1] - opt_point[d]);
            if (dist > max_dist) max_dist = dist;
        }
        
        if (max_dist > gkls_eps) isOk = false;

        if (isOk) {
            countSuccess++;
            totalIterations += (size_t)res[N + 1];
        } else {
            countFail++;
        }
    }

    // Останавливаем таймер
    auto end_time = std::chrono::high_resolution_clock::now();
    
    // Вычисляем разницу в секундах (с дробной частью)
    std::chrono::duration<double> diff = end_time - start_time;
    double totalTime = diff.count();
    double avgTime = totalTime / max_tests;

    double successRate = (double)countSuccess / max_tests * 100.0;
    double avgIter = countSuccess > 0 ? (double)totalIterations / countSuccess : 0.0;

    return {max_tests, countSuccess, countFail, successRate, avgIter, totalTime, avgTime};
}

inline void runGKLSTests() 
{
    vector<double> dists   = { 0.9, 0.9, 0.66, 0.9, 0.66, 0.9, 0.66, 0.66, 0.9, 0.9, 0.9 };
    vector<double> radiuss = { 0.2, 0.1, 0.2, 0.2, 0.2, 0.2, 0.3, 0.2, 0.2, 0.2, 0.2 };
    vector<int>    dims    = { 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 8 };
    vector<double> epss    = { 0.01, 0.01, 0.01, 0.01, 0.03, 0.03, 0.03, 0.03, 0.0464, 0.0517, 0.0562 };
    vector<string> shs     = { "simple", "hard", "simple", "hard", "simple", "hard", "simple", "hard", "simple", "simple", "simple" };
    vector<double> rs      = { 4.5, 5.6, 4.5, 5.6, 4.5, 5.6, 4.5, 5.6, 5.6, 5.6, 5.6 };
std::string directoryName = "results";
    
    std::filesystem::create_directories(directoryName);

    string filename = directoryName + "/GKLS_Results.csv";

    ofstream file(filename);
    if (file.is_open()) {
        file.imbue(std::locale(file.getloc(), new CommaNumPunct()));
        file << "Dim;Class;eps;r;Success_Rate_%;Success_Count;Fail_Count;Avg_Iterations;Total_Time_s;Avg_Time_s\n";
    }

    cout << "============================================================\n";
    cout << "  STARTING GKLS EXPERIMENTS (800 functions)\n";
    cout << "============================================================\n";

    for (int i = 0; i <= 7; i++) 
    {
        int N = dims[i];
        double gkls_eps = epss[i];
        string sh0 = shs[i];
        double r = rs[i];
        
        GKLSClass type = (sh0 == "simple") ? Simple : Hard;

        // int tightness = std::ceil(-std::log2(gkls_eps / 4.0));
        // if (tightness * N > 53) tightness = 53 / N; 
        int tightness = 7;
        cout << "GKLS | Dim=" << N << " | " << setw(6) << sh0 
             << " | r=" << fixed << setprecision(1) << r 
             << " | tight=" << tightness << " ... ";
        
        TestStats stats;
        
        if (N == 2) stats = run_gkls_batch<2>(gkls_eps, r, type, 100, tightness);
        else if (N == 3) stats = run_gkls_batch<3>(gkls_eps, r, type, 100, tightness);
        else if (N == 4) stats = run_gkls_batch<4>(gkls_eps, r, type, 100, tightness);
        else if (N == 5) stats = run_gkls_batch<5>(gkls_eps, r, type, 100, tightness);

        if (file.is_open()) {
            file << N << ";" << sh0 << ";" << gkls_eps << ";" << r << ";"
                 << stats.successRate << ";" << stats.success << ";"
                 << stats.fail << ";" << stats.avgIterations << ";"
                 << stats.totalTimeSec << ";" << stats.avgTimeSec << "\n";
        }
        
        cout << (stats.successRate == 100 ? "PERFECT!" : "Done.") 
             << " Rate: " << setw(5) << stats.successRate << "% | Iters: " 
             << fixed << setprecision(0) << stats.avgIterations 
             << " | Time: " << setprecision(2) << stats.totalTimeSec << "s\n";
    }
    
    cout << "============================================================\n";
    cout << "Эксперименты GKLS завершены. Результаты в " << filename << "\n\n";
    if (file.is_open()) file.close();
}

inline void runAllTests() 
{
    // eps для эвольвенты обычно берется более грубым, так как она растягивает отрезки
    double eps = 0.0001; 
    
    // Запускаем 100 функций Гришагина на отрезке [0, 1] по обеим координатам
    run_experiment_2d<TGrishaginProblem>(0.0, 1.0, eps, 100, "GRISHAGIN Problem", "Grishagin_Results.csv");
}