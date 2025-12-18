#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include "Algorithm.h"
#include <HillProblem.hpp>
#include <HillProblem.cpp>
using namespace std;

int main() 
{
  TPoint<double, 1> lowerBound(0);
  TPoint<double, 1> upperBound(1);
  double eps = 0.001;
  double r = 0;
  vector<double> res;
  for (size_t i = 0; i < 50; i++)
  {
    cout << i + 1 <<"\n";
    r = lConstantHill[i];
    THillProblem func(i);
    TAlgorithm<double, 1> a(lowerBound, upperBound, eps, r, &func);
    res = a.AGPStronginaMin();
    cout << minHill[i][0] - res[0] << '\n';
    cout << minHill[i][1] - res[1] << "\n";
  }
  return 0;
}