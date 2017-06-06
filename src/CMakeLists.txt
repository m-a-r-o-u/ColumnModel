add_library(columnmodel columnmodel.cpp)
add_executable(column  main.cpp)
target_link_libraries(column columnmodel)

