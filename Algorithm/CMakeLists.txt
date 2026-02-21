add_library(Algorithm INTERFACE)

target_include_directories(Algorithm INTERFACE 
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(Algorithm INTERFACE Point Interval OptProblem)