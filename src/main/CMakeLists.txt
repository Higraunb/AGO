add_executable(main main.cpp)

target_include_directories(main PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ../core/Algorithm
)

target_link_libraries(main PRIVATE Algorithm TestProblems OpenMP::OpenMP_CXX fmt::fmt spdlog::spdlog)