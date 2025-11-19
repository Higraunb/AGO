cmake_minimum_required(VERSION 3.16)

project(VisualAGP VERSION 0.1 LANGUAGES CXX)


if(MINGW)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static -static-libgcc -static-libstdc++")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
endif()

# Настройки Qt
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)


set(QWT_ROOT "C:/qwt-6.3.0")
set(QWT_BUILD_DIR "C:/qwt-6.3.0/build/Qt_5_12_12_mingw73_64-Debug")

# Прямое указание путей - закомментируйте то, что не подходит
if(EXISTS "${QWT_BUILD_DIR}/lib/libqwtd.a")
    # Для MinGW Debug
    set(QWT_LIBRARY "${QWT_BUILD_DIR}/lib/libqwtd.a")
elseif(EXISTS "${QWT_BUILD_DIR}/lib/libqwt.a") 
    # Для MinGW Release
    set(QWT_LIBRARY "${QWT_BUILD_DIR}/lib/libqwt.a")
elseif(EXISTS "${QWT_BUILD_DIR}/lib/qwtd.lib")
    # Для MSVC Debug
    set(QWT_LIBRARY "${QWT_BUILD_DIR}/lib/qwtd.lib")
elseif(EXISTS "${QWT_BUILD_DIR}/lib/qwt.lib")
    # Для MSVC Release  
    set(QWT_LIBRARY "${QWT_BUILD_DIR}/lib/qwt.lib")
else()
    # Поиск любого файла библиотеки
    file(GLOB QWT_LIB_FILES 
        "${QWT_BUILD_DIR}/lib/*qwt*.a"
        "${QWT_BUILD_DIR}/lib/*qwt*.lib"
        "${QWT_BUILD_DIR}/lib/*qwt*.dll.a"
    )
    if(QWT_LIB_FILES)
        list(GET QWT_LIB_FILES 0 QWT_LIBRARY)
    else()
        message(FATAL_ERROR "No Qwt library files found in ${QWT_BUILD_DIR}/lib/")
    endif()
endif()

set(QWT_INCLUDE_DIR "${QWT_ROOT}/src")

message(STATUS "QWT_LIBRARY: ${QWT_LIBRARY}")
message(STATUS "QWT_INCLUDE_DIR: ${QWT_INCLUDE_DIR}")

if(NOT EXISTS "${QWT_LIBRARY}")
    message(FATAL_ERROR "Qwt library file not found: ${QWT_LIBRARY}")
endif()

if(NOT EXISTS "${QWT_INCLUDE_DIR}/qwt_plot.h")
    message(FATAL_ERROR "Qwt headers not found in: ${QWT_INCLUDE_DIR}")
endif()

# Добавляем проект AGO как поддиректорию
add_subdirectory(AGO)  # Предполагая, что AGO находится в поддиректории с этим именем

# Получаем имена библиотек из AGO проекта
set(PROJECT_NAME AGO)
set(talgorithmlib TAlgorithmLibrary)
set(tfunctionlib TFunctionLibrary)
set(tintervallib TIntervalLibrary)
set(tpointlib TPointLibrary)

set(PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        mainwindow.h
        mainwindow.ui
)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(VisualAGP
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
    )
else()
    if(ANDROID)
        add_library(VisualAGP SHARED
            ${PROJECT_SOURCES}
        )
    else()
        add_executable(VisualAGP
            ${PROJECT_SOURCES}

        )
    endif()
endif()

# Подключаем библиотеки AGO к проекту
target_link_libraries(VisualAGP PRIVATE
    ${talgorithmlib}
    ${QWT_LIBRARY} 
    Qt${QT_VERSION_MAJOR}::Widgets
)

# Добавляем include директории из AGO
target_include_directories(VisualAGP PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/AGO/Algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/AGO/Function
    ${CMAKE_CURRENT_SOURCE_DIR}/AGO/Interval
    ${CMAKE_CURRENT_SOURCE_DIR}/AGO/Point
    ${QWT_INCLUDE_DIR}
)

target_compile_definitions(VisualAGP PRIVATE
    QWT_STATIC
)

# Остальные настройки...
set_target_properties(VisualAGP PROPERTIES
    ${BUNDLE_ID_OPTION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

include(GNUInstallDirs)
install(TARGETS VisualAGP
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(VisualAGP)
endif()
