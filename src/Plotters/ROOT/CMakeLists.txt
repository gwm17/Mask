find_package(ROOT REQUIRED)
add_executable(RootPlot)
target_include_directories(RootPlot 
    PUBLIC ${MASK_INCLUDE_DIR}
    SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS}
)

target_sources(RootPlot PUBLIC
    RootPlotter.cpp
)

target_link_libraries(RootPlot
    Mask
    ${ROOT_LIBRARIES}
)

set_target_properties(RootPlot PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${MASK_BINARY_DIR})