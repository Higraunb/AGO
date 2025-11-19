 file(GLOB hdrs "*.h*") #Добавляем в переменную hdrs все файлы с расширением .h
file(GLOB srcs "*.cpp")#Добавляем в переменную srcs все файлы с расширением .h

add_library(${talgorithmlib} STATIC ${srcs} ${hdrs})
target_include_directories(${talgorithmlib} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${talgorithmlib} ${tpointlib} ${tintervallib} ${tfunctionlib})