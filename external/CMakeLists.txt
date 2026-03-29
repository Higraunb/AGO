
file(GLOB OPT_PROBLEM_SRC 
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

add_library(OptProblem ${OPT_PROBLEM_SRC})


target_include_directories(OptProblem PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)


file(GLOB_RECURSE TEST_PROBLEMS_SRC 
    "${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.h"
)

add_library(TestProblems ${TEST_PROBLEMS_SRC})


target_include_directories(TestProblems PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Hill
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Shekel
)

target_link_libraries(TestProblems PUBLIC OptProblem)