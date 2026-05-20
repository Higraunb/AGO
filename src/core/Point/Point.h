#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <array>
#include "../../Logger/Logger.h"

template<std::size_t N>
class TPoint
{
private:
  std::array<double, N> coords;

public:
  TPoint();
  explicit TPoint(const std::vector<double>& values);
  explicit TPoint(std::vector<double>&& values);
  explicit TPoint(const double (&arr)[N]);
  template<class... Args>
  explicit TPoint(Args... args);

  double& operator[](std::size_t i);
  const double& operator[](std::size_t i) const;

  constexpr std::size_t size() const;

  std::vector<double> toVector() const;

  TPoint operator+(const TPoint& other) const;
  TPoint operator-(const TPoint& other) const;
  TPoint operator*(double scalar) const;
  TPoint operator/(double scalar) const;
  TPoint& operator+=(const TPoint& other);
  TPoint& operator-=(const TPoint& other);
  TPoint& operator*=(double scalar);
  TPoint& operator/=(double scalar);

  TPoint operator-() const;

  bool operator==(const TPoint& other) const;
  bool operator!=(const TPoint& other) const;
  bool operator<(const TPoint& other) const;

  double dot(const TPoint& other) const;

  template<std::size_t M = N>
  typename std::enable_if<M == 3, TPoint>::type cross(const TPoint& other) const;

  double norm() const;
  double distance(const TPoint& other) const;

  void resize(std::size_t newSize) { coords.resize(newSize); }

  TPoint normalized() const;
  void normalize();

  const std::vector<double>& data() const;
  std::vector<double> data();

  template<std::size_t M>
  friend std::ostream& operator << (std::ostream& out, const TPoint<M>& other);
};


template<std::size_t N>
TPoint<N>::TPoint() 
{
  coords.fill(0);
}

template<std::size_t N>
std::vector<double> TPoint<N>::toVector() const {
  return std::vector<double>(coords.begin(), coords.end());
}

template<std::size_t N>
TPoint<N>::TPoint(const double (&arr)[N])
{
  std::copy(std::begin(arr), std::end(arr), std::begin(coords));
}

template<std::size_t N>
TPoint<N>::TPoint(const std::vector<double>& values) 
{
  if (values.size() != N) 
    throw std::invalid_argument("Vector size does not match point dimension");
  std::copy(values.begin(), values.end(), coords.begin());
}

template<std::size_t N>
TPoint<N>::TPoint(std::vector<double>&& values) : coords(std::move(values)) 
{}

template<std::size_t N>
template<class... Args>
TPoint<N>::TPoint(Args... args) : coords{ static_cast<double>(args)... } 
{
  static_assert(sizeof...(args) == N, "Wrong number of arguments");
}

template<std::size_t N>
double& TPoint<N>::operator[](std::size_t i)
{
  if (i >= N) {
      LOG_ERROR("TPoint::operator[] - Выход за границы! Индекс: {}, Макс: {}", i, N - 1);
      throw std::out_of_range("Index out of bounds");
  }
  return coords[i];
}

template<std::size_t N>
const double& TPoint<N>::operator[](std::size_t i) const
{
  if (i >= N) {
      LOG_ERROR("TPoint::operator[] const - Выход за границы! Индекс: {}, Макс: {}", i, N - 1);
      throw std::out_of_range("Index out of bounds");
  }
  return coords[i];
}

template<std::size_t N>
constexpr std::size_t TPoint<N>::size() const 
{ 
  return N; 
}

template<std::size_t N>
TPoint<N> TPoint<N>::operator+(const TPoint& other) const 
{
  TPoint result;
  for (std::size_t i = 0; i < N; ++i)
    result[i] = coords[i] + other[i];

  return result;
}

template<std::size_t N>
TPoint<N> TPoint<N>::operator-(const TPoint& other) const 
{
  TPoint result;
  for (std::size_t i = 0; i < N; ++i)
    result[i] = coords[i] - other[i];

  return result;
}

template<std::size_t N>
TPoint<N> TPoint<N>::operator*(double scalar) const 
{
  TPoint result;
  for (std::size_t i = 0; i < N; ++i)
    result[i] = coords[i] * scalar;

  return result;
}

template<std::size_t N>
TPoint<N> TPoint<N>::operator/(double scalar) const 
{
  if (std::abs(scalar) < 1e-10)
    throw std::domain_error("Division by zero");
  TPoint result;
  for (std::size_t i = 0; i < N; ++i)
    result[i] = coords[i] / scalar;

  return result;
}

template<std::size_t N>
TPoint<N>& TPoint<N>::operator+=(const TPoint& other) 
{
  for (std::size_t i = 0; i < N; ++i)
    coords[i] += other[i];

  return *this;
}

template<std::size_t N>
TPoint<N>& TPoint<N>::operator-=(const TPoint& other) 
{
  for (std::size_t i = 0; i < N; ++i)
    coords[i] -= other[i];

  return *this;
}

template<std::size_t N>
TPoint<N>& TPoint<N>::operator*=(double scalar)
{
  for (std::size_t i = 0; i < N; ++i) 
    coords[i] *= scalar;

  return *this;
}

template<std::size_t N>
TPoint<N>& TPoint<N>::operator/=(double scalar) 
{
  if (std::abs(scalar) < 1e-10) 
    throw std::domain_error("Division by zero");

  for (std::size_t i = 0; i < N; ++i)
    coords[i] /= scalar;

  return *this;
}

template<std::size_t N>
TPoint<N> TPoint<N>::operator-() const 
{
  TPoint result;
  for (std::size_t i = 0; i < N; ++i) 
    result[i] = -coords[i];

  return result;
}

template<std::size_t N>
bool TPoint<N>::operator==(const TPoint& other) const 
{
  for (std::size_t i = 0; i < N; ++i) 
    if (std::abs(coords[i] - other[i]) > 1e-10) 
      return false;
  return true;
}

template<std::size_t N>
bool TPoint<N>::operator!=(const TPoint& other) const 
{
  return !(*this == other);
}

template<std::size_t N>
bool TPoint<N>::operator<(const TPoint& other) const
{
  return coords < other.coords;
}


template<std::size_t N>
double TPoint<N>::dot(const TPoint& other) const {
  double result = 0;
  for (std::size_t i = 0; i < N; ++i)
    result += coords[i] * other[i];

  return result;
}

template<std::size_t N>
double TPoint<N>::norm() const { return std::sqrt(dot(*this)); }

template<std::size_t N>
template<std::size_t M>
typename std::enable_if<M == 3, TPoint<N>>::type
TPoint<N>::cross(const TPoint& other) const 
{
  return TPoint(
    coords[1] * other[2] - coords[2] * other[1],
    coords[2] * other[0] - coords[0] * other[2],
    coords[0] * other[1] - coords[1] * other[0]
  );
}

template<std::size_t N>
double TPoint<N>::distance(const TPoint& other) const { return (*this - other).norm(); }

template<std::size_t N>
TPoint<N> TPoint<N>::normalized() const
{
  double n = norm();
  if (n < 1e-10) 
  {
    LOG_ERROR("TPoint::normalized - Попытка нормализовать нулевой вектор!");
    throw std::domain_error("Cannot normalize zero vector");
  }
  return *this / n;
}

template<std::size_t N>
void TPoint<N>::normalize() 
{
  *this = normalized();
}

template<std::size_t N>
const std::vector<double>& TPoint<N>::data() const 
{ 
  return coords; 
}

template<std::size_t N>
std::vector<double> TPoint<N>::data() 
{ 
  return coords;
}

template<std::size_t N>
TPoint<N> operator*(double scalar, const TPoint<N>& p) 
{
  return p * scalar;
}

template<std::size_t N>
std::ostream& operator<<(std::ostream& out, const TPoint<N>& other)
{
  out << "(";
  for (std::size_t i = 0; i < N; ++i) 
  {
    out << other.coords[i];
    if (i < N - 1)
      out << ", ";
  }
  out << ")";
  return out;
}
