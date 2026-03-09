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
    T zLeft;
    T zRight;
    T vLeft;
    T vRight;

    Interval(T l, T r, T zl, T zr, T vl, T vr) : left(l), right(r), zLeft(zl), zRight(zr), R(0), vRight(vr), vLeft(vl) {}
    
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

  void initialize(T left, T right, T zLeft, T zRight, T vLeft = 0, T vRight = 0);

  std::size_t splitByIndex(std::size_t index, T x, T zx, T vx);

  const std::vector<Interval>& getIntervals() const;
  std::vector<Interval>& getIntervals();

  typename std::vector<Interval>::iterator findInterval(T x);
  typename std::vector<Interval>::const_iterator findInterval(T x) const;

  Interval getMaxRInterval() const;
  T getRight(size_t index) const{return intervals[index].right;};
  T getLeft(size_t index) const{return intervals[index].left;};
  T getLength(size_t index) const{return intervals[index].length();};
  T getZLeft(size_t index) const { return intervals[index].zLeft; }
  T getZRight(size_t index) const { return intervals[index].zRight; }
  T getVLeft(size_t index) const { return intervals[index].vLeft; }
  T getVRight(size_t index) const { return intervals[index].vRight; }
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
  initialize(0, 1, 1, 1);
}

template<class T>
TInterval<T>::TInterval(T left, T right) 
{
  initialize(left, right);
}

template<class T>
void TInterval<T>::initialize(T left, T right, T zLeft, T zRight, T vLeft, T vRight) 
{
  if (left >= right) {
    throw std::invalid_argument("Left boundary must be less than right boundary");
  }
  intervals.clear();
  intervals.push_back(Interval(left, right, zLeft, zRight, vLeft, vRight));
}

template<class T>
std::size_t TInterval<T>::splitByIndex(std::size_t index, T x, T zx, T vx) 
  {
    if (index >= intervals.size()) throw std::out_of_range("Index out of range");

    Interval& old = intervals[index];
    T right = old.right;
    T zRight = old.zRight;
    T vRight = old.vRight;

    old.right = x;
    old.zRight = zx;
    old.vRight = vx;

    intervals.push_back(Interval(x, right, zx, zRight, vx, vRight));

    

    return intervals.size() - 1; 
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