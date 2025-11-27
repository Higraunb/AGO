#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include "Algorithm.h"
using namespace std;

double testFunction(const TPoint<double, 1>& point) 
{
  double x = point[0];
  double res = (x * x);
  return res;
}
int main() 
{
  TPoint<double, 1> lowerBound(-1.8);
  TPoint<double, 1> upperBound(2.2);
  double eps = 0.01;
  double r = 2;
  // TAlgorithm<double, 1> a(lowerBound, upperBound, eps, r, testFunction);
  // a.AGPStronginaMax();
  return 0;
}