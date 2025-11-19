#include <array>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <sstream>

template<class T, std::size_t N>
class TPoint {
private:
  std::array<T, N> coords;

public:
  TPoint();
  explicit TPoint(const std::array<T, N>& values);
  explicit TPoint(std::array<T, N>&& values);

  template<class... Args>
  explicit TPoint(Args... args);

  T& operator[](std::size_t i);
  const T& operator[](std::size_t i) const;

  constexpr std::size_t size() const;

  auto begin();
  auto end();
  auto begin() const;
  auto end() const;

  TPoint operator+(const TPoint& other) const;
  TPoint operator-(const TPoint& other) const;
  TPoint operator*(T scalar) const;
  TPoint operator/(T scalar) const;
  TPoint& operator+=(const TPoint& other);
  TPoint& operator-=(const TPoint& other);
  TPoint& operator*=(T scalar);
  TPoint& operator/=(T scalar);

  TPoint operator-() const;

  bool operator==(const TPoint& other) const;
  bool operator!=(const TPoint& other) const;
  bool operator<(const TPoint& other) const;

  T dot(const TPoint& other) const;

  template<std::size_t M = N>
  typename std::enable_if<M == 3, TPoint>::type cross(const TPoint& other) const;

  T squaredNorm() const;
  T norm() const;
  T distance(const TPoint& other) const;

  TPoint normalized() const;
  void normalize();

  TPoint projectOn(const TPoint& other) const;

  const std::array<T, N>& data() const;
  std::array<T, N>& data();

  T to_T();
  template<class O, std::size_t M>
  friend std::ostream& operator << (std::ostream& out, const TPoint<O, M>& other);
};


template<class T, std::size_t N>
TPoint<T, N>::TPoint() : coords{}
{
}

template<class T, std::size_t N>
TPoint<T, N>::TPoint(const std::array<T, N>& values) : coords(values)
{
}

template<class T, std::size_t N>
TPoint<T, N>::TPoint(std::array<T, N>&& values) : coords(std::move(values))
{
}

template<class T, std::size_t N>
template<class... Args>
TPoint<T, N>::TPoint(Args... args) : coords{ static_cast<T>(args)... }
{
  static_assert(sizeof...(args) == N, "Wrong number of arguments");
}

template<class T, std::size_t N>
T& TPoint<T, N>::operator[](std::size_t i)
{
  if (i >= N) throw
    std::out_of_range("Index out of range");
  return coords[i];
}

template<class T, std::size_t N>
const T& TPoint<T, N>::operator[](std::size_t i) const
{
  if (i >= N) throw
    std::out_of_range("Index out of range");
  return coords[i];
}

template<class T, std::size_t N>
constexpr std::size_t TPoint<T, N>::size() const
{
  return N;
}

template<class T, std::size_t N>
auto TPoint<T, N>::begin()
{
  return coords.begin();
}

template<class T, std::size_t N>
auto TPoint<T, N>::end()
{
  return coords.end();
}

template<class T, std::size_t N>
auto TPoint<T, N>::begin() const
{
  return coords.begin();
}

template<class T, std::size_t N>
auto TPoint<T, N>::end() const
{
  return coords.end();
}

template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::operator+(const TPoint& other) const
{
  TPoint result;
  for (std::size_t i = 0; i < N; ++i)
    result[i] = coords[i] + other[i];

  return result;
}

template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::operator-(const TPoint& other) const
{
  TPoint result;
  for (std::size_t i = 0; i < N; ++i)
    result[i] = coords[i] - other[i];

  return result;
}

template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::operator*(T scalar) const
{
  TPoint result;
  for (std::size_t i = 0; i < N; ++i)
    result[i] = coords[i] * scalar;

  return result;
}

template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::operator/(T scalar) const {
  if (std::abs(scalar) < 1e-10) {
    throw std::domain_error("Division by zero");
  }
  TPoint result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = coords[i] / scalar;
  }
  return result;
}

template<class T, std::size_t N>
TPoint<T, N>& TPoint<T, N>::operator+=(const TPoint& other)
{
  for (std::size_t i = 0; i < N; ++i)
    coords[i] += other[i];

  return *this;
}

template<class T, std::size_t N>
TPoint<T, N>& TPoint<T, N>::operator-=(const TPoint& other)
{
  for (std::size_t i = 0; i < N; ++i)
    coords[i] -= other[i];

  return *this;
}

template<class T, std::size_t N>
TPoint<T, N>& TPoint<T, N>::operator*=(T scalar)
{
  for (std::size_t i = 0; i < N; ++i)
    coords[i] *= scalar;

  return *this;
}

template<class T, std::size_t N>
TPoint<T, N>& TPoint<T, N>::operator/=(T scalar)
{
  if (std::abs(scalar) < 1e-10)
    throw std::domain_error("Division by zero");

  for (std::size_t i = 0; i < N; ++i)
    coords[i] /= scalar;

  return *this;
}

template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::operator-() const
{
  TPoint result;
  for (std::size_t i = 0; i < N; ++i)
    result[i] = -coords[i];

  return result;
}

template<class T, std::size_t N>
bool TPoint<T, N>::operator==(const TPoint& other) const
{
  for (std::size_t i = 0; i < N; ++i) {
    if (std::abs(coords[i] - other[i]) > 1e-10) {
      return false;
    }
  }
  return true;
}

template<class T, std::size_t N>
bool TPoint<T, N>::operator!=(const TPoint& other) const
{
  return !(*this == other);
}

template<class T, std::size_t N>
bool TPoint<T, N>::operator<(const TPoint& other) const
{
  return coords < other.coords;
}

template<class T, std::size_t N>
T TPoint<T, N>::squaredNorm() const
{
    T result = 0;
    for (std::size_t i = 0; i < N; ++i) {
        result += coords[i] * coords[i];
    }
    return result;
}

template<class T, std::size_t N>
T TPoint<T, N>::dot(const TPoint& other) const {
  T result = 0;
  for (std::size_t i = 0; i < N; ++i) {
    result += coords[i] * other[i];
  }
  return result;
}

template<class T, std::size_t N>
template<std::size_t M>
typename std::enable_if<M == 3, TPoint<T, N>>::type
TPoint<T, N>::cross(const TPoint& other) const {
  return TPoint(
    coords[1] * other[2] - coords[2] * other[1],
    coords[2] * other[0] - coords[0] * other[2],
    coords[0] * other[1] - coords[1] * other[0]
  );
}

template<class T, std::size_t N>
T TPoint<T, N>::norm() const
{
  return std::sqrt(squaredNorm());
}

template<class T, std::size_t N>
T TPoint<T, N>::distance(const TPoint& other) const {
  return (*this - other).norm();
}

template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::normalized() const
{
  T n = norm();
  if (n < 1e-10)
    throw std::domain_error("Cannot normalize zero vector");
  return *this / n;
}

template<class T, std::size_t N>
void TPoint<T, N>::normalize()
{
  *this = normalized();
}

template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::projectOn(const TPoint& other) const
{
  T scale = this->dot(other) / other.squaredNorm();
  return other * scale;
}

template<class T, std::size_t N>
const std::array<T, N>& TPoint<T, N>::data() const
{
  return coords;
}

template<class T, std::size_t N>
std::array<T, N>& TPoint<T, N>::data()
{
  return coords;
}

template<class T, std::size_t N>
inline T TPoint<T, N>::to_T()
{
  return coords[0];
}

template<class T, std::size_t N>
TPoint<T, N> operator*(T scalar, const TPoint<T, N>& p)
{
  return p * scalar;
}

template<class O, std::size_t M>
std::ostream& operator<<(std::ostream& out, const TPoint<O, M>& other)
{
  out << "(";
  for (std::size_t i = 0; i < M; ++i)
  {
    out << other.coords[i];
    if (i < M - 1)
      out << ", ";
  }
  out << ")";
  return out;
}
