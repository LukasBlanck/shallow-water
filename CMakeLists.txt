cmake_minimum_required(VERSION 3.16)
project(my_project CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(my_project
    apps/main.cpp
    src/io/netcdf_writer.cpp
    configs/simulation_config.cpp
    configs/config.hpp
)

target_include_directories(my_project
    PRIVATE
    ${CMAKE_SOURCE_DIR}
)

find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET netcdf)
pkg_check_modules(NETCDFCXX REQUIRED IMPORTED_TARGET netcdf-cxx4)

target_link_libraries(my_project
    PRIVATE
    PkgConfig::NETCDF
    PkgConfig::NETCDFCXX
)