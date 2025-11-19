 file(GLOB hdrs "*.h*") #Добавляем в переменную hdrs все файлы с расширением .h
file(GLOB srcs "*.cpp")#Добавляем в переменную srcs все файлы с расширением .h

add_library(${tintervallib} STATIC ${srcs} ${hdrs})
target_include_directories(${tintervallib} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${tintervallib} ${tpointlib})