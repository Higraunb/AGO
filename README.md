# AGO Application (GSA)

![C++](https://img.shields.io/badge/C++-17-blue.svg) ![CMake](https://img.shields.io/badge/CMake-3.14+-green.svg)

Консольное приложение для запуска бенчмарков и тестирования эффективности алгоритмов глобального поиска (включая функции GKLS). 

## Структура проекта

Проект разделен на несколько основных модулей:
* `src/core/` — базовые примитивы и логика (Point, Interval, Algorithm).
* `src/main/` — точка входа (`main.cpp`) и модуль запуска тестов (`TestRunner`).
* `Logger/` — интерфейс логирования.
* `external/` — сторонние библиотеки и тестовые классы задач.

## Зависимости

Для сборки проекта требуются:
* **C++17**
* **CMake 3.14** или новее
* **OpenCV**
* **Threads**

Следующие библиотеки скачиваются и собираются автоматически через `FetchContent`:
* **fmt** (10.2.1)
* **spdlog** (v1.13.0)

## Сборка (Linux / WSL)
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
