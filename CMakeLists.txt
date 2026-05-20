cmake_minimum_required(VERSION 3.14)
project(AGO CXX)

# ==============================================================================
# 1. Принудительная установка Release сборки
# ==============================================================================
# Записываем "Release" в кэш CMake принудительно (FORCE), чтобы гарантировать
# срабатывание флагов оптимизации при любой генерации проекта.
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)

# Задаем стандарт C++17 (можно поменять на 20, если нужно)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ==============================================================================
# 2. Флаги компилятора
# ==============================================================================
add_compile_options(-march=native -ffast-math)

# ==============================================================================
# 3. Подключение сторонних библиотек (FetchContent)
# ==============================================================================
include(FetchContent)

# Скачиваем fmt
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG        10.2.1
)
FetchContent_MakeAvailable(fmt)

# Скачиваем spdlog
FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG        v1.13.0
)
FetchContent_MakeAvailable(spdlog)

# Ищем OpenCV в системе (собирать его из исходников каждый раз слишком долго)
set(OpenCV_DIR "C:/msys64/ucrt64/lib/cmake/opencv4")
find_package(OpenCV REQUIRED)

# Ищем системную библиотеку потоков (нужно для Linux, на Windows пройдет незаметно)
find_package(Threads REQUIRED)

# ==============================================================================
# 4. Сборка проекта
# ==============================================================================

# Создаем библиотеку Logger (header-only, так как только Logger.h)
add_library(AGO_Logger INTERFACE)
target_include_directories(AGO_Logger INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Logger)
target_link_libraries(AGO_Logger INTERFACE spdlog::spdlog)

# Подключаем все поддиректории
add_subdirectory(external)
add_subdirectory(src/core/Point)
add_subdirectory(src/core/Interval)
add_subdirectory(src/core/Algorithm)
add_subdirectory(src/main)