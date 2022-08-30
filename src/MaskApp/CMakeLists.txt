add_executable(MaskApp)
target_include_directories(MaskApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_sources(MaskApp PUBLIC
    main.cpp
)

target_link_libraries(MaskApp
    Mask
)

set_target_properties(MaskApp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${MASK_BINARY_DIR})