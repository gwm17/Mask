add_executable(Kinematics)
target_include_directories(Kinematics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_sources(Kinematics PUBLIC
    main.cpp
)

target_link_libraries(Kinematics
    Mask
)

set_target_properties(Kinematics PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${MASK_BINARY_DIR})