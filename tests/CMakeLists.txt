# ==========================================
# ЦЕЛЬ 1: Библиотека интерфейсов (OptProblem)
# Собираем файлы из папки src/
# ==========================================
file(GLOB OPT_PROBLEM_SRC 
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

add_library(OptProblem ${OPT_PROBLEM_SRC})

# Делаем папку src доступной для всех, кто подключит эту либу
target_include_directories(OptProblem PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)


# ==========================================
# ЦЕЛЬ 2: Библиотека задач (TestProblems)
# Собираем файлы из папки sample_src/
# ==========================================
file(GLOB_RECURSE TEST_PROBLEMS_SRC 
    "${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.h"
)

add_library(TestProblems ${TEST_PROBLEMS_SRC})

# Указываем пути к заголовкам
target_include_directories(TestProblems PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Hill
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Shekel
)

# Задачи (Hill, Shekel) наследуются от IOptProblem, поэтому линкуем первую цель
target_link_libraries(TestProblems PUBLIC OptProblem)