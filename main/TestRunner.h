#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <locale>
#include <cmath>
#include "Algorithm.h"

// Включаем функции Гришагина
#include "grishagin_function.hpp"

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
        vector<double> res = alg.Solve(10000, true); 

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
    if (countFail > 0 && countFail <= 3)
    {
        cout << "  [!] Внимание: при r = " << fixed << setprecision(1) << r << " упало тестов: " << countFail << "\n";
        for (const auto& f : failures)
        {
            cout << "      Тест " << setw(3) << f.id 
                 << " | Ожидали: (" << setw(5) << f.expX << ", " << setw(5) << f.expY << ")"
                 << " | Получили: (" << setw(5) << f.gotX << ", " << setw(5) << f.gotY << ")"
                 << " | Z_Ожид: " << setw(8) << f.expZ << " | Z_Получ: " << setw(8) << f.gotZ << "\n";
        }
    }

    double successRate = (double)countSuccess / max_tests * 100.0;
    double avgIter = countSuccess > 0 ? (double)totalIterations / countSuccess : 0.0;

    return {max_tests, countSuccess, countFail, successRate, avgIter};
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
    int tightness = 10; 

    for (double r = 1.5; r <= 8.1; r += 0.2) 
    {
        int tightness = 3;
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

inline void runAllTests() 
{
    // eps для эвольвенты обычно берется более грубым, так как она растягивает отрезки
    double eps = 0.0001; 
    
    // Запускаем 100 функций Гришагина на отрезке [0, 1] по обеим координатам
    run_experiment_2d<TGrishaginProblem>(0.0, 1.0, eps, 100, "GRISHAGIN Problem", "Grishagin_Results.csv");
}