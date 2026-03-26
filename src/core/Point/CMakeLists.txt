add_library(Point
    Point.h
    Point.cpp
)

target_include_directories(Point PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Убедитесь, что Point доступен для других целей
set(Point_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)