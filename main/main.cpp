#include "Algorithm.h"

double func(TPoint<double, 1> other)
{
  double x = other[0];
  return x * x;
}
int main()
{
  TPoint<double, 1> low(-2.0), higt(2.0);
  TAlgorithm<double, 1> a(low, higt, 0.001, 2, func);
  a.AGPStronginaMax();
  return 0;
}