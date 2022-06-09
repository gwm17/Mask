cmake_minimum_required(VERSION 3.16)

project(Mask)

set(MASK_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(MASK_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(MASK_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src/vendor/catima)
add_subdirectory(src/Mask)
add_subdirectory(src/MaskApp)
add_subdirectory(src/Detectors)
add_subdirectory(src/Plotters/ROOT)