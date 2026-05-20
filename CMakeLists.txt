cmake_minimum_required(VERSION 3.14)
project(AGO CXX)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-march=native -ffast-math)

include(FetchContent)

FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG        10.2.1
)
FetchContent_MakeAvailable(fmt)

FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG        v1.13.0
)
FetchContent_MakeAvailable(spdlog)

set(OpenCV_DIR "C:/msys64/ucrt64/lib/cmake/opencv4")
find_package(OpenCV REQUIRED)

find_package(Threads REQUIRED)

add_library(AGO_Logger INTERFACE)
target_include_directories(AGO_Logger INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Logger)
target_link_libraries(AGO_Logger INTERFACE spdlog::spdlog)

add_subdirectory(external)
add_subdirectory(src/core/Point)
add_subdirectory(src/core/Interval)
add_subdirectory(src/core/Algorithm)
add_subdirectory(src/main)