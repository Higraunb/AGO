add_library(Algorithm
    Algorithm.h
    Algorithm.cpp
)

target_include_directories(Algorithm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Algorithm PUBLIC Point Interval Function)