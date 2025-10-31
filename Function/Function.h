#include <functional>
#include <stdexcept>
#include <cmath>
#include "Interval.h"

template<class T, std::size_t N>
class TFunction {
public:
  using FunctionType = std::function<T(const TPoint<T, N>&)>;

private:
  FunctionType func;
  TPoint<T, N> lowerBound;
  TPoint<T, N> upperBound;

public:
  TFunction();
  TFunction(FunctionType f, const TPoint<T, N>& lower, const TPoint<T, N>& upper);

  // Установить функцию и границы
  void setFunction(FunctionType f);
  void setBounds(const TPoint<T, N>& lower, const TPoint<T, N>& upper);

  // Нормировка: из реальных координат [lower, upper] в [0, 1]^N
  TPoint<T, N> normalize(const TPoint<T, N>& point) const;

  // Обратная нормировка: из [0, 1]^N в реальные координаты [lower, upper]
  TPoint<T, N> denormalize(const TPoint<T, N>& point) const;

  // Вычислить значение функции в точке (в реальных координатах)
  T operator()(const TPoint<T, N>& point) const;

  // Вычислить значение функции в нормированной точке [0, 1]^N
  T evaluateNormalized(const TPoint<T, N>& normalizedPoint) const;

  // Получить границы
  const TPoint<T, N>& getLowerBound() const;
  const TPoint<T, N>& getUpperBound() const;

  // Размерность
  constexpr std::size_t dimension() const;
};


template<class T, std::size_t N>
TFunction<T, N>::TFunction() 
  : func(nullptr)
{}

template<class T, std::size_t N>
TFunction<T, N>::TFunction(FunctionType f, const TPoint<T, N>& lower, const TPoint<T, N>& upper)
  : func(f), lowerBound(lower), upperBound(upper)
{
  // Проверяем корректность границ
  for (std::size_t i = 0; i < N; ++i) {
    if (lowerBound[i] >= upperBound[i]) {
      throw std::invalid_argument("Lower bound must be less than upper bound in all dimensions");
    }
  }
}

template<class T, std::size_t N>
void TFunction<T, N>::setFunction(FunctionType f)
{
  func = f;
}

template<class T, std::size_t N>
void TFunction<T, N>::setBounds(const TPoint<T, N>& lower, const TPoint<T, N>& upper)
{
  for (std::size_t i = 0; i < N; ++i) {
    if (lower[i] >= upper[i]) {
      throw std::invalid_argument("Lower bound must be less than upper bound in all dimensions");
    }
  }
  lowerBound = lower;
  upperBound = upper;
}

template<class T, std::size_t N>
TPoint<T, N> TFunction<T, N>::normalize(const TPoint<T, N>& point) const
{
  TPoint<T, N> normalized;
  for (std::size_t i = 0; i < N; ++i) {
    T range = upperBound[i] - lowerBound[i];
    if (std::abs(range) < 1e-10) {
      throw std::domain_error("Zero range in dimension");
    }
    normalized[i] = (point[i] - lowerBound[i]) / range;
  }
  return normalized;
}

template<class T, std::size_t N>
TPoint<T, N> TFunction<T, N>::denormalize(const TPoint<T, N>& point) const
{
  TPoint<T, N> denormalized;
  for (std::size_t i = 0; i < N; ++i) {
    T range = upperBound[i] - lowerBound[i];
    denormalized[i] = lowerBound[i] + point[i] * range;
  }
  return denormalized;
}

template<class T, std::size_t N>
T TFunction<T, N>::operator()(const TPoint<T, N>& point) const
{
  if (!func) {
    throw std::runtime_error("Function is not set");
  }
  return func(point);
}

template<class T, std::size_t N>
T TFunction<T, N>::evaluateNormalized(const TPoint<T, N>& normalizedPoint) const
{
  if (!func) {
    throw std::runtime_error("Function is not set");
  }
  
  // Преобразуем из [0, 1]^N в реальные координаты и вычисляем
  TPoint<T, N> realPoint = denormalize(normalizedPoint);
  return func(realPoint);
}

template<class T, std::size_t N>
const TPoint<T, N>& TFunction<T, N>::getLowerBound() const
{
  return lowerBound;
}

template<class T, std::size_t N>
const TPoint<T, N>& TFunction<T, N>::getUpperBound() const
{
  return upperBound;
}

template<class T, std::size_t N>
constexpr std::size_t TFunction<T, N>::dimension() const
{
  return N;
}