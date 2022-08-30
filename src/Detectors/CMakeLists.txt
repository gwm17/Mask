add_executable(Detectors)
target_include_directories(Detectors PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_sources(Detectors PUBLIC
    AnasenDeadChannelMap.cpp
    AnasenDeadChannelMap.h
    AnasenArray.cpp
    AnasenArray.h
    main.cpp
    DetectorArray.h
    QQQDetector.cpp
    QQQDetector.h
    SabreDeadChannelMap.cpp
    SabreDeadChannelMap.h
    SabreDetector.cpp
    SabreDetector.h
    SabreArray.cpp
    SabreArray.h
    SX3Detector.cpp
    SX3Detector.h
)

target_link_libraries(Detectors
    Mask
)

set_target_properties(Detectors PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${MASK_BINARY_DIR})