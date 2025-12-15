add_library(TFunction
    Function.h
    Function.cpp
)

target_include_directories(TFunction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TFunction PUBLIC Point Interval)