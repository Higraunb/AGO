add_library(Function
    Function.h
    Function.cpp
)

target_include_directories(Function PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Function PUBLIC Point Interval)