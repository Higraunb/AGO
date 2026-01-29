#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>

template<class T>
class TInterval 
{
public:
  struct Interval 
  {
    T left;
    T right;
    T R;  
    
    Interval(T l, T r) : left(l), right(r), R(0) {}
    
    T length() const { return right - left; }
    T center() const { return (left + right) / 2; }
    
    bool operator<(const Interval& other) const {
      return left < other.left;
    }
  };

private:
  std::vector<Interval> intervals;

public:
  TInterval();
  explicit TInterval(T left, T right);

  void initialize(T left, T right);

  std::pair<Interval, Interval> split(T x);

  const std::vector<Interval>& getIntervals() const;
  std::vector<Interval>& getIntervals();

  typename std::vector<Interval>::iterator findInterval(T x);
  typename std::vector<Interval>::const_iterator findInterval(T x) const;

  Interval getMaxRInterval() const;
  T getRight(size_t index) const;
  T getLeft(size_t index) const;
  T getLength(size_t index) const;
  std::size_t getMaxRIntervalIndex() const;
  std::size_t getMinRIntervalIndex() const;

  void setIntervalR(size_t index, T R);

  std::size_t size() const;

  void clear();
  
  bool empty() const;

  template<class U>
  friend std::ostream& operator<<(std::ostream& out, const TInterval<U>& ti);
};


template<class T>
TInterval<T>::TInterval() 
{
  initialize(0, 1);
}

template<class T>
TInterval<T>::TInterval(T left, T right) 
{
  initialize(left, right);
}

template<class T>
void TInterval<T>::initialize(T left, T right) 
{
  if (left >= right) {
    throw std::invalid_argument("Left boundary must be less than right boundary");
  }
  intervals.clear();
  intervals.push_back(Interval(left, right));
}

template<class T>
std::pair<typename TInterval<T>::Interval, typename TInterval<T>::Interval> 
TInterval<T>::split(T x) 
{
  auto it = findInterval(x);
  
  if (it == intervals.end()) {
    throw std::out_of_range("Point x is not within any interval");
  }

  T left = it->left;
  T right = it->right;

  if (std::abs(x - left) < 1e-10 || std::abs(x - right) < 1e-10) {
    throw std::invalid_argument("Point x coincides with interval boundary");
  }

  Interval left_interval(left, x);
  Interval right_interval(x, right);

  std::size_t idx = std::distance(intervals.begin(), it);

  intervals.erase(it);
  intervals.insert(intervals.begin() + idx, {left_interval, right_interval});

  return {left_interval, right_interval};
}

template<class T>
inline const std::vector<typename TInterval<T>::Interval>& 
TInterval<T>::getIntervals() const 
{
  return intervals;
}

template<class T>
inline std::vector<typename TInterval<T>::Interval>& 
TInterval<T>::getIntervals()
{
  return intervals;
}

template<class T>
typename std::vector<typename TInterval<T>::Interval>::iterator 
inline TInterval<T>::findInterval(T x) 
{
  for (auto it = intervals.begin(); it != intervals.end(); ++it) {
    if (x >= it->left && x <= it->right) {
      return it;
    }
  }
  
  return intervals.end();
}

template<class T>
typename std::vector<typename TInterval<T>::Interval>::const_iterator 
inline TInterval<T>::findInterval(T x) const 
{
  for (auto it = intervals.begin(); it != intervals.end(); ++it) {
    if (x >= it->left && x <= it->right) {
      return it;
    }
  }
  return intervals.end();
}

template<class T>
typename TInterval<T>::Interval 
inline TInterval<T>::getMaxRInterval() const 
{
  if (intervals.empty()) {
    throw std::runtime_error("No intervals available");
  }

  return *std::max_element(intervals.begin(), intervals.end(),
    [](const Interval& a, const Interval& b) {
      return a.R < b.R;
    });
}

template <class T>
inline T TInterval<T>::getRight(size_t index) const
{
  return intervals[index].right;
}

template <class T>
inline T TInterval<T>::getLeft(size_t index) const
{
  return intervals[index].left;
}

template <class T>
inline T TInterval<T>::getLength(size_t index) const
{
  return intervals[index].length();
}

template<class T>
inline std::size_t TInterval<T>::getMaxRIntervalIndex() const 
{
  if (intervals.empty()) {
    throw std::runtime_error("No intervals available");
  }

  auto it = std::max_element(intervals.begin(), intervals.end(),
    [](const Interval& a, const Interval& b) {
      return a.R < b.R;
    });

  return std::distance(intervals.begin(), it);
}

template<class T>
inline std::size_t TInterval<T>::getMinRIntervalIndex() const 
{
  if (intervals.empty())
    throw std::runtime_error("No intervals available");

  auto it = std::min_element(intervals.begin(), intervals.end(),
  [](const Interval& a, const Interval& b) {
    return a.R < b.R;
  });

  return std::distance(intervals.begin(), it);
}

template <class T>
inline void TInterval<T>::setIntervalR(size_t index, T R)
{
  if(index >= intervals.size())
    throw std::out_of_range("Index out of range");
  intervals[index].R = R;
}

template<class T>
inline std::size_t TInterval<T>::size() const 
{
  return intervals.size();
}

template<class T>
inline void TInterval<T>::clear() 
{
  intervals.clear();
}

template<class T>
inline bool TInterval<T>::empty() const 
{
  return intervals.empty();
}

template<class T>
std::ostream& operator<<(std::ostream& out, const TInterval<T>& ti) 
{
  out << "TInterval with " << ti.size() << " intervals:\n";
  for (const auto& interval : ti.intervals) {
    out << "  [" << interval.left << ", " << interval.right 
        << "] length=" << interval.length() 
        << " R=" << interval.R << "\n";
  }
  return out;
}