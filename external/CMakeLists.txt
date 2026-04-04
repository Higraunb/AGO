# Собираем все .cpp исходники из папки external рекурсивно
file(GLOB_RECURSE EXTERNAL_SOURCES
    "src/*.cpp"
    "sample_src/*.cpp"
)

# Создаем статическую библиотеку
add_library(AGO_External STATIC ${EXTERNAL_SOURCES})

# Указываем пути для поиска заголовочных файлов.
# PUBLIC означает, что тот, кто слинкует AGO_External, тоже получит эти пути.
target_include_directories(AGO_External PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Hill
    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Shekel
)