add_library(Algorithm Algorithm.cpp)

target_include_directories(Algorithm PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(Algorithm PUBLIC Point Interval OptProblem OpenMP::OpenMP_CXX)