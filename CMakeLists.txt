cmake_minimum_required(VERSION 3.16)


if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE "Release")
    message("Building release")
else()
    message("Building debug")
endif()

project(Mask)

set(MASK_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(MASK_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(MASK_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED COMPONENTS GenVector)

add_subdirectory(src/vendor/catima)
set(YAML_CPP_BUILD_TOOLS Off CACHE BOOL "Enable parse tools")
add_subdirectory(src/vendor/yaml-cpp)
add_subdirectory(src/Mask)
add_subdirectory(src/Kinematics)
add_subdirectory(src/Detectors)
add_subdirectory(src/Plotters)