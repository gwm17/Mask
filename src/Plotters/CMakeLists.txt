add_executable(RootPlot)
target_include_directories(RootPlot 
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
           ${CMAKE_CURRENT_SOURCE_DIR}/..
    SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS}
)

target_sources(RootPlot PUBLIC
    RootPlotter.cpp
    RootPlotter.h
    main.cpp
)

target_link_libraries(RootPlot
    Mask
    ${ROOT_LIBRARIES}
)

set_target_properties(RootPlot PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${MASK_BINARY_DIR})