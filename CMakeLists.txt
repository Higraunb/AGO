cmake_minimum_required(VERSION 3.10)
project(AGO)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_subdirectory(Point)
add_subdirectory(Interval)

# Сначала собираем тесты, так как там лежат интерфейсы (IOptProblem)
add_subdirectory(tests)

# Потом алгоритм, который их использует
add_subdirectory(Algorithm)

add_subdirectory(main)