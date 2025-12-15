# Папки src и sample_src находятся ВНУТРИ tests, поэтому используем CMAKE_CURRENT_SOURCE_DIR
file(GLOB SRC_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Hill/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Shekel/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp
)

# Отладка - выведем найденные файлы
message(STATUS "Source files for tests library:")
foreach(file ${SRC_FILES})
    message(STATUS "  ${file}")
endforeach()

# Проверка что файлы найдены
if(NOT SRC_FILES)
    message(FATAL_ERROR "No source files found in ${CMAKE_CURRENT_SOURCE_DIR}/sample_src or ${CMAKE_CURRENT_SOURCE_DIR}/src")
endif()

add_library(tests ${SRC_FILES})

# Публичные include директории для библиотеки
target_include_directories(tests PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Hill
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Shekel
)

# Линкуем с другими библиотеками проекта
target_link_libraries(tests PUBLIC 
    Algorithm
    TFunction
    Interval
    Point
)