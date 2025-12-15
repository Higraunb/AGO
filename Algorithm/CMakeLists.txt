add_library(Algorithm
    Algorithm.h
    Algorithm.cpp
)

target_include_directories(Algorithm PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/sample_src
)

target_link_libraries(Algorithm PUBLIC Point Interval TFunction tests)