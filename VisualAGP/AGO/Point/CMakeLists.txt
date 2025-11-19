 file(GLOB hdrs "*.h*") #Добавляем в переменную hdrs все файлы с расширением .h
file(GLOB srcs "*.cpp")#Добавляем в переменную srcs все файлы с расширением .h

add_library(${tpointlib} STATIC ${srcs} ${hdrs})
target_include_directories(${tpointlib} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})