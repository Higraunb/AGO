cmake_minimum_required(VERSION 3.10)
project(AGO)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Ищем OpenMP в системе
find_package(OpenMP REQUIRED)

# Ищем fmt и spdlog в системе
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)

# Внешние библиотеки (оптимизационные задачи)
add_subdirectory(external)

# Основные компоненты проекта
add_subdirectory(src/core/Point)
add_subdirectory(src/core/Interval)
add_subdirectory(src/core/Algorithm)
add_subdirectory(src/main)