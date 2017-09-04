find_package(yaml-cpp)
include_directories(${YAML_CPP_INCLUDE_DIR})

find_package(fpda_rrtm)
include_directories(${FPDA_RRTM_INCLUDE_DIR})

find_package(netCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIR})

add_library(columnmodel 
            thermodynamic.cpp
            tau_relax.cpp
            columnmodel.cpp)

        target_link_libraries(columnmodel ${YAML_CPP_LIBRARIES} ${FPDA_RRTM_LIBRARIES} ${NETCDF_LIBRARIES} netcdf_c++4)

target_include_directories(columnmodel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/)

add_executable(column
               main.cpp
               )
target_link_libraries(column columnmodel)
