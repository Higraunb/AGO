file(GLOB_RECURSE EXTERNAL_SOURCES
    "src/*.cpp"
    "sample_src/*.cpp"
)

add_library(AGO_External STATIC ${EXTERNAL_SOURCES})

target_include_directories(AGO_External PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Hill
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Shekel
)

target_link_libraries(AGO_External PUBLIC AGO_Logger)
