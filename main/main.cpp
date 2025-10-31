#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <Function.h>
using namespace std;

double testFunction(const TPoint<double, 1>& point) 
{
  double x = point[0];
  double res = (sin(x) * pow(2.71828182, x) - 7);
  return res;
}
int main() 
{
  TPoint<double, 1> lowerBound({-2.0});
  TPoint<double, 1> upperBound({2.0});
  TFunction<double, 1> func(testFunction, lowerBound, upperBound);
  TInterval<double> interval;
  double eps;
  double L = 0;
  eps = 0.01;
  double r;
  r = 2;
  double R = 0;
  do
  {
    double mu = 0;
    for (size_t i = 0; i < interval.size(); i++)
    {
      double zi = func.evaluateNormalized(TPoint<double,1>(interval.getRigth(i)));
      double zi1 = func.evaluateNormalized(TPoint<double,1>(interval.getLeft(i)));
      double xi = interval.getRigth(i);
      double xi1 = interval.getLeft(i);
      mu = max(mu, abs(zi - zi1) / (xi - xi1));
    }
    if (mu == 0)
      L = 1;
    else
      L = r * mu;
    for (size_t i = 0; i < interval.size(); i++)
    {
      double zi = func.evaluateNormalized(TPoint<double,1>(interval.getRigth(i)));
      double zi1 = func.evaluateNormalized(TPoint<double,1>(interval.getLeft(i)));
      double xi = interval.getRigth(i);
      double xi1 = interval.getLeft(i);
      R = L * (xi - xi1) + (pow(zi - zi1, 2) / (L * (xi - xi1))) - 2 * (zi + zi1);
      interval.setIntervalR(i, R);
    }
    size_t maxIndex = interval.getMaxRIntervalIndex();
    double zi = func.evaluateNormalized(TPoint<double,1>(interval.getRigth(maxIndex)));
    double zi1 = func.evaluateNormalized(TPoint<double,1>(interval.getLeft(maxIndex)));
    double xi = interval.getRigth(maxIndex);
    double xi1 = interval.getLeft(maxIndex);
    double newX = ((xi1 + xi) / 2) - ((zi - zi1) / (2 * L));
    interval.split(newX);
  } while (interval.getLength(interval.getMaxRIntervalIndex()) > eps);
  size_t maxIndex = interval.getMaxRIntervalIndex();
  double x0 = interval.getLeft(maxIndex);
  double x1 = interval.getRigth(maxIndex);
  TPoint<double, 1> xr = (func.denormalize(TPoint<double, 1>(x0)) + func.denormalize(TPoint<double, 1>(x1)))/2;
  cout << "xr :" << xr << endl;
  cout << "func :" << func(xr) << endl; 
  return 0;
}