cmake_minimum_required(VERSION 3.15)
project(my_project CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_OPENMP "Build with OpenMP support if available" ON)
option(REQUIRE_OPENMP "Fail configuration if OpenMP is not available" OFF)

add_executable(my_project
    apps/main.cpp
    src/io/netCDF_writer.cpp
    configs/simulation_config.cpp
    tests/sanity_checks/sanity_checks.cpp
    include/io/sanity_checks_netcdf_writer.cpp
)

target_include_directories(my_project
    PRIVATE
    ${CMAKE_SOURCE_DIR}
)

find_package(PkgConfig QUIET)

if(PkgConfig_FOUND)
    pkg_check_modules(NETCDF QUIET IMPORTED_TARGET netcdf)
    pkg_check_modules(NETCDFCXX QUIET IMPORTED_TARGET netcdf-cxx4)
endif()

if(NETCDF_FOUND)
    message(STATUS "Found NetCDF C via pkg-config")
    target_link_libraries(my_project PRIVATE PkgConfig::NETCDF)
else()
    find_path(NETCDF_INCLUDE_DIR_FOUND
        NAMES netcdf.h
        PATHS
            ${NETCDF_INCLUDE_DIR}
            $ENV{NETCDF_INCLUDE_DIR}
            $ENV{UIBK_NETCDF_4_INC}
        NO_DEFAULT_PATH
    )

    find_library(NETCDF_LIBRARY_FOUND
        NAMES netcdf
        PATHS
            ${NETCDF_LIBRARY_DIR}
            $ENV{NETCDF_LIBRARY_DIR}
            $ENV{UIBK_NETCDF_4_LIB}
        NO_DEFAULT_PATH
    )

    if(NOT NETCDF_INCLUDE_DIR_FOUND OR NOT NETCDF_LIBRARY_FOUND)
        message(FATAL_ERROR "Could not find NetCDF C. Need netcdf.h and libnetcdf.")
    endif()

    message(STATUS "Found NetCDF C include: ${NETCDF_INCLUDE_DIR_FOUND}")
    message(STATUS "Found NetCDF C library: ${NETCDF_LIBRARY_FOUND}")

    target_include_directories(my_project PRIVATE ${NETCDF_INCLUDE_DIR_FOUND})
    target_link_libraries(my_project PRIVATE ${NETCDF_LIBRARY_FOUND})
endif()

if(NETCDFCXX_FOUND)
    message(STATUS "Found NetCDF C++ via pkg-config")
    target_link_libraries(my_project PRIVATE PkgConfig::NETCDFCXX)
else()
    find_path(NETCDFCXX_INCLUDE_DIR_FOUND
        NAMES ncFile.h netcdf
        PATHS
            ${NETCDF_CXX_INCLUDE_DIR}
            $ENV{NETCDF_CXX_INCLUDE_DIR}
        NO_DEFAULT_PATH
    )

    find_library(NETCDFCXX_LIBRARY_FOUND
        NAMES netcdf_c++4 netcdf-cxx4
        PATHS
            ${NETCDF_CXX_LIBRARY_DIR}
            $ENV{NETCDF_CXX_LIBRARY_DIR}
        NO_DEFAULT_PATH
    )

    if(NOT NETCDFCXX_INCLUDE_DIR_FOUND OR NOT NETCDFCXX_LIBRARY_FOUND)
        message(FATAL_ERROR "Could not find NetCDF C++. Need ncFile.h and libnetcdf_c++4.")
    endif()

    message(STATUS "Found NetCDF C++ include: ${NETCDFCXX_INCLUDE_DIR_FOUND}")
    message(STATUS "Found NetCDF C++ library: ${NETCDFCXX_LIBRARY_FOUND}")

    target_include_directories(my_project PRIVATE ${NETCDFCXX_INCLUDE_DIR_FOUND})
    target_link_libraries(my_project PRIVATE ${NETCDFCXX_LIBRARY_FOUND})
endif()

# ============================================================
# Optional OpenMP
# ============================================================

if(ENABLE_OPENMP)
    list(APPEND CMAKE_PREFIX_PATH
        "$ENV{OpenMP_ROOT}"
        "$ENV{OPENMP_ROOT}"
        "$ENV{HOMEBREW_PREFIX}"
        "/opt/homebrew/opt/libomp"
        "/usr/local/opt/libomp"
        "/opt/homebrew/opt/llvm"
        "/usr/local/opt/llvm"
    )

    find_package(OpenMP QUIET COMPONENTS CXX)

    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found")
        message(STATUS "OpenMP_CXX_FLAGS=${OpenMP_CXX_FLAGS}")
        target_link_libraries(my_project PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(my_project PRIVATE USE_OPENMP=1)
    elseif(REQUIRE_OPENMP)
        message(FATAL_ERROR
            "OpenMP was requested but not found.\n"
            "On macOS/Homebrew try:\n"
            "  brew install libomp\n"
            "  cmake -S . -B build -DOpenMP_ROOT=/opt/homebrew/opt/libomp\n"
            "or:\n"
            "  brew install llvm libomp\n"
            "  cmake -S . -B build "
            "-DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang "
            "-DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++\n"
            "To build without OpenMP:\n"
            "  cmake -S . -B build -DENABLE_OPENMP=OFF"
        )
    else()
        message(WARNING
            "OpenMP not found. Building serial/non-OpenMP solver only.\n"
            "To require OpenMP, configure with -DREQUIRE_OPENMP=ON."
        )
        target_compile_definitions(my_project PRIVATE USE_OPENMP=0)
    endif()
else()
    message(STATUS "OpenMP disabled")
    target_compile_definitions(my_project PRIVATE USE_OPENMP=0)
endif()