#pragma once

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "../../Logger/Logger.h"

class Intervals 
{
public:
  struct Interval 
  {
    double left;
    double right;
    double R;  
    double zLeft;
    double zRight;
    double vLeft;
    double vRight;

    Interval(double l, double r, double zl, double zr, double vl, double vr) 
      : left(l), right(r), zLeft(zl), zRight(zr), R(0.0), vLeft(vl), vRight(vr) {}
    
    double length() const { return right - left; }
    double center() const { return (left + right) / 2.0; }
    
    bool operator<(const Interval& other) const { return left < other.left; }
  };

private:
  std::vector<Interval> intervals;

public:
  Intervals();
  explicit Intervals(double left, double right);

  void initialize(double left, double right, double zLeft, double zRight, double vLeft = 0.0, double vRight = 0.0);

  std::size_t splitByIndex(std::size_t index, double x, double zx, double vx);

  const std::vector<Interval>& getIntervals() const;
  std::vector<Interval>& getIntervals();

  std::vector<Interval>::iterator findInterval(double x);
  std::vector<Interval>::const_iterator findInterval(double x) const;

  Interval getMaxRInterval() const;
  double getRight(size_t index) const { return intervals[index].right; }
  double getLeft(size_t index) const { return intervals[index].left; }
  double getLength(size_t index) const { return intervals[index].length(); }
  double getZLeft(size_t index) const { return intervals[index].zLeft; }
  double getZRight(size_t index) const { return intervals[index].zRight; }
  double getVLeft(size_t index) const { return intervals[index].vLeft; }
  double getVRight(size_t index) const { return intervals[index].vRight; }
  std::size_t getMaxRIntervalIndex() const;
  std::size_t getMinRIntervalIndex() const;

  void setIntervalR(size_t index, double R);

  std::size_t size() const;

  void clear();
  
  bool empty() const;

  friend std::ostream& operator<<(std::ostream& out, const Intervals& ti);
};

inline Intervals::Intervals() 
{
  initialize(0.0, 1.0, 1.0, 1.0);
}

inline Intervals::Intervals(double left, double right) 
{
  initialize(left, right, 1.0, 1.0);
}

inline void Intervals::initialize(double left, double right, double zLeft, double zRight, double vLeft, double vRight) 
{
  if (left >= right)
    throw std::invalid_argument("Left boundary must be less than right boundary");
    
  intervals.clear();
  intervals.push_back(Interval(left, right, zLeft, zRight, vLeft, vRight));
}

inline std::size_t Intervals::splitByIndex(std::size_t index, double x, double zx, double vx) 
{
  if (index >= intervals.size()) throw std::out_of_range("Index out of range");

  Interval& old = intervals[index];
  double right = old.right;
  double zRight = old.zRight;
  double vRight = old.vRight;

  old.right = x;
  old.zRight = zx;
  old.vRight = vx;

  intervals.push_back(Interval(x, right, zx, zRight, vx, vRight));

  return intervals.size() - 1; 
}

inline const std::vector<Intervals::Interval>& Intervals::getIntervals() const 
{
  return intervals;
}

inline std::vector<Intervals::Interval>& Intervals::getIntervals()
{
  return intervals;
}

inline std::vector<Intervals::Interval>::iterator Intervals::findInterval(double x) 
{
  for (auto it = intervals.begin(); it != intervals.end(); ++it) {
    if (x >= it->left && x <= it->right) {
      return it;
    }
  }
  return intervals.end();
}

inline std::vector<Intervals::Interval>::const_iterator Intervals::findInterval(double x) const 
{
  for (auto it = intervals.begin(); it != intervals.end(); ++it) {
    if (x >= it->left && x <= it->right) {
      return it;
    }
  }
  return intervals.end();
}

inline Intervals::Interval Intervals::getMaxRInterval() const 
{
  if (intervals.empty()) {
    LOG_ERROR("Intervals::getMaxRInterval - Нет доступных интервалов (вектор пуст)!");
    throw std::runtime_error("No intervals available");
  }

  return *std::max_element(intervals.begin(), intervals.end(),
    [](const Interval& a, const Interval& b) {
      return a.R < b.R;
    });
}

inline std::size_t Intervals::getMaxRIntervalIndex() const 
{
  if(intervals.empty()) {
    LOG_ERROR("Intervals::getMaxRIntervalIndex - Нет доступных интервалов (вектор пуст)!");
    throw std::out_of_range("No intervals available");
  }

  auto it = std::max_element(intervals.begin(), intervals.end(),
    [](const Interval& a, const Interval& b) {
      return a.R < b.R;
    });

  return std::distance(intervals.begin(), it);
}

inline std::size_t Intervals::getMinRIntervalIndex() const 
{
  if(intervals.empty())
  {
    LOG_ERROR("Intervals::getMinRIntervalIndex - Нет доступных интервалов (вектор пуст)!");
    throw std::out_of_range("No intervals available");
  }

  auto it = std::min_element(intervals.begin(), intervals.end(),
  [](const Interval& a, const Interval& b) {
    return a.R < b.R;
  });

  return std::distance(intervals.begin(), it);
}

inline void Intervals::setIntervalR(size_t index, double R)
{
  if(index >= intervals.size())
    throw std::out_of_range("Index out of range");
  intervals[index].R = R;
}

inline std::size_t Intervals::size() const 
{
  return intervals.size();
}

inline void Intervals::clear() 
{
  intervals.clear();
}

inline bool Intervals::empty() const 
{
  return intervals.empty();
}

inline std::ostream& operator<<(std::ostream& out, const Intervals& ti) 
{
  out << "Intervals with " << ti.size() << " intervals:\n";
  for (const auto& interval : ti.intervals) 
  {
    out << "  [" << interval.left << ", " << interval.right 
        << "] length=" << interval.length() 
        << " R=" << interval.R << "\n";
  }
  return out;
}