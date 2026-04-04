# Создаем исполняемый файл
add_executable(AGO_Main main.cpp)

# Добавляем папку main в инклуды (чтобы main.cpp видел TestRunner.h), а также пути OpenCV
target_include_directories(AGO_Main PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)

# Линкуем все наши модули и новые библиотеки
target_link_libraries(AGO_Main PRIVATE
    AGO_External
    AGO_Point
    AGO_Interval
    AGO_Algorithm
    AGO_Logger
    Threads::Threads
    fmt::fmt
    spdlog::spdlog
    ${OpenCV_LIBS}
)