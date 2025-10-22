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
  // Конструкторы
  TPoint();
  explicit TPoint(const std::array<T, N>& values);
  explicit TPoint(std::array<T, N>&& values);

  template<class... Args>
  explicit TPoint(Args... args);

  // Доступ к элементам
  T& operator[](std::size_t i);
  const T& operator[](std::size_t i) const;
  T& at(std::size_t i);
  const T& at(std::size_t i) const;

  // Размерность
  constexpr std::size_t dimension() const;
  constexpr std::size_t size() const;

  // Итераторы
  auto begin();
  auto end();
  auto begin() const;
  auto end() const;

  // Арифметические операции
  TPoint operator+(const TPoint& other) const;
  TPoint operator-(const TPoint& other) const;
  TPoint operator*(T scalar) const;
  TPoint operator/(T scalar) const;
  TPoint& operator+=(const TPoint& other);
  TPoint& operator-=(const TPoint& other);
  TPoint& operator*=(T scalar);
  TPoint& operator/=(T scalar);

  // Унарный минус
  TPoint operator-() const;

  // Операции сравнения
  bool operator==(const TPoint& other) const;
  bool operator!=(const TPoint& other) const;
  bool operator<(const TPoint& other) const;

  // Векторные операции
  T dot(const TPoint& other) const;

  template<std::size_t M = N>
  typename std::enable_if<M == 3, TPoint>::type cross(const TPoint& other) const;

  // Нормы и расстояния
  T norm() const;
  T squaredNorm() const;
  T manhattanNorm() const;
  T maxNorm() const;
  T distance(const TPoint& other) const;
  T squaredDistance(const TPoint& other) const;

  // Нормализация
  TPoint normalized() const;
  void normalize();

  // Проекция
  TPoint projectOn(const TPoint& other) const;

  // Строковое представление
  std::string toString() const;

  // Доступ к данным
  const std::array<T, N>& data() const;
  std::array<T, N>& data();
};

// Определения конструкторов
template<class T, std::size_t N>
TPoint<T, N>::TPoint() : coords{} {}

template<class T, std::size_t N>
TPoint<T, N>::TPoint(const std::array<T, N>& values) : coords(values) {}

template<class T, std::size_t N>
TPoint<T, N>::TPoint(std::array<T, N>&& values) : coords(std::move(values)) {}

template<class T, std::size_t N>
template<class... Args>
TPoint<T, N>::TPoint(Args... args) : coords{ static_cast<T>(args)... } {
  static_assert(sizeof...(args) == N, "Wrong number of arguments");
}

// Определения методов доступа
template<class T, std::size_t N>
T& TPoint<T, N>::operator[](std::size_t i) {
  if (i >= N) throw std::out_of_range("Index out of range");
  return coords[i];
}

template<class T, std::size_t N>
const T& TPoint<T, N>::operator[](std::size_t i) const {
  if (i >= N) throw std::out_of_range("Index out of range");
  return coords[i];
}

template<class T, std::size_t N>
T& TPoint<T, N>::at(std::size_t i) {
  return coords.at(i);
}

template<class T, std::size_t N>
const T& TPoint<T, N>::at(std::size_t i) const {
  return coords.at(i);
}

// Определения методов размерности
template<class T, std::size_t N>
constexpr std::size_t TPoint<T, N>::dimension() const { return N; }

template<class T, std::size_t N>
constexpr std::size_t TPoint<T, N>::size() const { return N; }

// Определения итераторов
template<class T, std::size_t N>
auto TPoint<T, N>::begin() { return coords.begin(); }

template<class T, std::size_t N>
auto TPoint<T, N>::end() { return coords.end(); }

template<class T, std::size_t N>
auto TPoint<T, N>::begin() const { return coords.begin(); }

template<class T, std::size_t N>
auto TPoint<T, N>::end() const { return coords.end(); }

// Определения арифметических операций
template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::operator+(const TPoint& other) const {
  TPoint result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = coords[i] + other[i];
  }
  return result;
}

template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::operator-(const TPoint& other) const {
  TPoint result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = coords[i] - other[i];
  }
  return result;
}

template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::operator*(T scalar) const {
  TPoint result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = coords[i] * scalar;
  }
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
TPoint<T, N>& TPoint<T, N>::operator+=(const TPoint& other) {
  for (std::size_t i = 0; i < N; ++i) {
    coords[i] += other[i];
  }
  return *this;
}

template<class T, std::size_t N>
TPoint<T, N>& TPoint<T, N>::operator-=(const TPoint& other) {
  for (std::size_t i = 0; i < N; ++i) {
    coords[i] -= other[i];
  }
  return *this;
}

template<class T, std::size_t N>
TPoint<T, N>& TPoint<T, N>::operator*=(T scalar) {
  for (std::size_t i = 0; i < N; ++i) {
    coords[i] *= scalar;
  }
  return *this;
}

template<class T, std::size_t N>
TPoint<T, N>& TPoint<T, N>::operator/=(T scalar) {
  if (std::abs(scalar) < 1e-10) {
    throw std::domain_error("Division by zero");
  }
  for (std::size_t i = 0; i < N; ++i) {
    coords[i] /= scalar;
  }
  return *this;
}

// Определение унарного минуса
template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::operator-() const {
  TPoint result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = -coords[i];
  }
  return result;
}

// Определения операций сравнения
template<class T, std::size_t N>
bool TPoint<T, N>::operator==(const TPoint& other) const {
  for (std::size_t i = 0; i < N; ++i) {
    if (std::abs(coords[i] - other[i]) > 1e-10) {
      return false;
    }
  }
  return true;
}

template<class T, std::size_t N>
bool TPoint<T, N>::operator!=(const TPoint& other) const {
  return !(*this == other);
}

template<class T, std::size_t N>
bool TPoint<T, N>::operator<(const TPoint& other) const {
  return coords < other.coords;
}

// Определения векторных операций
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

// Определения норм и расстояний
template<class T, std::size_t N>
T TPoint<T, N>::norm() const {
  return std::sqrt(squaredNorm());
}

template<class T, std::size_t N>
T TPoint<T, N>::squaredNorm() const {
  T sum = 0;
  for (const auto& c : coords) {
    sum += c * c;
  }
  return sum;
}

template<class T, std::size_t N>
T TPoint<T, N>::manhattanNorm() const {
  T sum = 0;
  for (const auto& c : coords) {
    sum += std::abs(c);
  }
  return sum;
}

template<class T, std::size_t N>
T TPoint<T, N>::maxNorm() const {
  T maxVal = std::abs(coords[0]);
  for (std::size_t i = 1; i < N; ++i) {
    maxVal = std::max(maxVal, std::abs(coords[i]));
  }
  return maxVal;
}

template<class T, std::size_t N>
T TPoint<T, N>::distance(const TPoint& other) const {
  return (*this - other).norm();
}

template<class T, std::size_t N>
T TPoint<T, N>::squaredDistance(const TPoint& other) const {
  return (*this - other).squaredNorm();
}

// Определения нормализации
template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::normalized() const {
  T n = norm();
  if (n < 1e-10) {
    throw std::domain_error("Cannot normalize zero vector");
  }
  return *this / n;
}

template<class T, std::size_t N>
void TPoint<T, N>::normalize() {
  *this = normalized();
}

// Определение проекции
template<class T, std::size_t N>
TPoint<T, N> TPoint<T, N>::projectOn(const TPoint& other) const {
  T scale = this->dot(other) / other.squaredNorm();
  return other * scale;
}

// Определение строкового представления
template<class T, std::size_t N>
std::string TPoint<T, N>::toString() const {
  std::ostringstream oss;
  oss << "(";
  for (std::size_t i = 0; i < N; ++i) {
    oss << coords[i];
    if (i < N - 1) oss << ", ";
  }
  oss << ")";
  return oss.str();
}

// Определения доступа к данным
template<class T, std::size_t N>
const std::array<T, N>& TPoint<T, N>::data() const { return coords; }

template<class T, std::size_t N>
std::array<T, N>& TPoint<T, N>::data() { return coords; }

// Внешние операторы
template<class T, std::size_t N>
TPoint<T, N> operator*(T scalar, const TPoint<T, N>& p) {
  return p * scalar;
}

template<class T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const TPoint<T, N>& p) {
  os << p.toString();
  return os;
}

// Удобные псевдонимы
template<class T>
using Point1D = TPoint<T, 1>;

template<class T>
using Point2D = TPoint<T, 2>;

template<class T>
using Point3D = TPoint<T, 3>;

template<class T>
using Point4D = TPoint<T, 4>;

// Специализации для популярных типов
using Point2f = Point2D<float>;
using Point2d = Point2D<double>;
using Point3f = Point3D<float>;
using Point3d = Point3D<double>;
using Point4f = Point4D<float>;
using Point4d = Point4D<double>;