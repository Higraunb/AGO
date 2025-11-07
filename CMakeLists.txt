cmake_minimum_required(VERSION 3.10...3.31) #Проверка версии CMake.

set(PROJECT_NAME AGO) #Создаем обычную локальную переменную с именем проекта
project(${PROJECT_NAME}) # Название проекта


set(application AGOAPP) #Переменная с именем приложения

set(talgorithmlib TAlgorithmLibrary)
set(tfunctionlib TFunctionLibrary)
set(tintervallib TIntervalLibrary)
set(tpointlib TPointLibrary)
   

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)
if(NOT CMAKE_BUILD_TYPE) # проверка имеющихся у компилятора вариантов сборки
  set(CMAKE_BUILD_TYPE Release) #Если нет то устанавливаем Release
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Algorithm)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Function)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Interval)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Point)

add_subdirectory(Algorithm)
add_subdirectory(Function)
add_subdirectory(Interval)
add_subdirectory(Point)
add_subdirectory(main)
add_subdirectory(tests)
