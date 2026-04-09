cmake_minimum_required(VERSION 3.14)
project(AGO CXX)

# Задаем стандарт C++17 (можно поменять на 20, если нужно)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    # -O3 включает самые агрессивные оптимизации
    # -march=native оптимизирует под твой конкретный процессор
    # -ffast-math может чуть-чуть снизить точность флоатов в крайних случаях, но сильно ускоряет математику
    add_compile_options(-O3 -march=native -ffast-math)
endif()

# Подключаем модуль для автоматического скачивания библиотек
include(FetchContent)

# 1. Скачиваем fmt
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG        10.2.1
)
FetchContent_MakeAvailable(fmt)

# 2. Скачиваем spdlog
FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG        v1.13.0
)
FetchContent_MakeAvailable(spdlog)

# 3. Ищем OpenCV в системе (собирать его из исходников каждый раз слишком долго)
set(OpenCV_DIR "C:/msys64/ucrt64/lib/cmake/opencv4")
find_package(OpenCV REQUIRED)

# Ищем системную библиотеку потоков (нужно для Linux, на Windows пройдет незаметно)
find_package(Threads REQUIRED)

# 2. Создаем библиотеку Logger (header-only, так как только Logger.h)
add_library(AGO_Logger INTERFACE)
target_include_directories(AGO_Logger INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Logger)
target_link_libraries(AGO_Logger INTERFACE spdlog::spdlog)

# 3. Подключаем все поддиректории
add_subdirectory(external)
add_subdirectory(src/core/Point)
add_subdirectory(src/core/Interval)
add_subdirectory(src/core/Algorithm)
add_subdirectory(src/main)