add_executable(DetectEff)
target_include_directories(DetectEff PUBLIC ${MASK_INCLUDE_DIR})

target_sources(DetectEff PUBLIC
    AnasenDeadChannelMap.cpp
    AnasenEfficiency.cpp
    DetectorEfficiency.cpp
    QQQDetector.cpp
    SabreDeadChannelMap.cpp
    SabreDetector.cpp
    SabreEfficiency.cpp
    StripDetector.cpp
)

target_link_libraries(DetectEff
    Mask
)

set_target_properties(DetectEff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${MASK_BINARY_DIR})