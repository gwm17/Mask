add_executable(DetectEff)
target_include_directories(DetectEff PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_sources(DetectEff PUBLIC
    AnasenDeadChannelMap.cpp
    AnasenDeadChannelMap.h
    AnasenEfficiency.cpp
    AnasenEfficiency.h
    main.cpp
    DetectorEfficiency.h
    QQQDetector.cpp
    QQQDetector.h
    SabreDeadChannelMap.cpp
    SabreDeadChannelMap.h
    SabreDetector.cpp
    SabreDetector.h
    SabreEfficiency.cpp
    SabreEfficiency.h
    StripDetector.cpp
    StripDetector.h
)

target_link_libraries(DetectEff
    MaskDict
    Mask
    catima
)

set_target_properties(DetectEff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${MASK_BINARY_DIR})