find_package(yaml-cpp)
include_directories(${YAML_CPP_INCLUDE_DIR})

add_library(columnmodel 
            state.cpp
            thermodynamic.cpp
            columnmodel.cpp)

target_include_directories(columnmodel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/)

add_executable(column
               main.cpp
               )
target_link_libraries(column columnmodel)
