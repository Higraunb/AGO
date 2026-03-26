add_library(Interval
    Interval.h
    Interval.cpp
)

target_include_directories(Interval PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Убедитесь, что Interval доступен для других целей
set(Interval_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)