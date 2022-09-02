add_library(MaskDict SHARED)
target_include_directories(MaskDict
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS}
)

ROOT_GENERATE_DICTIONARY(mask_dict Nucleus.h LINKDEF LinkDef_Nucleus.h MODULE MaskDict)

target_sources(MaskDict PRIVATE Nucleus.h Nucleus.cpp MassLookup.h MassLookup.cpp)
target_link_libraries(MaskDict ${ROOT_LIBRARIES})
set_target_properties(MaskDict PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${MASK_LIBRARY_DIR})
add_custom_command(TARGET MaskDict POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy
	${CMAKE_CURRENT_BINARY_DIR}/libMaskDict_rdict.pcm
	${MASK_LIBRARY_DIR}/libMaskDict_rdict.pcm
	)

add_library(Mask STATIC)
target_include_directories(Mask 
    PUBLIC ${MASK_INCLUDE_DIR}
    ${ROOT_INCLUDE_DIRS}
)

target_sources(Mask PRIVATE
    AngularDistribution.cpp
    AngularDistribution.h
    DecaySystem.cpp
    DecaySystem.h
    LayeredTarget.cpp
    LayeredTarget.h
    LegendrePoly.cpp
    LegendrePoly.h
    MaskApp.cpp
    MaskApp.h
    MassLookup.cpp
    MassLookup.h
    OneStepSystem.cpp
    OneStepSystem.h
    RandomGenerator.cpp
    RandomGenerator.h
    Reaction.cpp
    Reaction.h
    ReactionSystem.cpp
    ReactionSystem.h
    Stopwatch.cpp
    Stopwatch.h
    Target.cpp
    Target.h
    ThreeStepSystem.cpp
    ThreeStepSystem.h
    TwoStepSystem.cpp
    TwoStepSystem.h
    ThreadPool.h
    FileWriter.h
    FileWriter.cpp
)

set(THREADS_PREFER_PTHREAD_FLAG On)
find_package(Threads REQUIRED)

target_link_libraries(Mask catima MaskDict ${ROOT_LIBRARIES} Threads::Threads)
set_target_properties(Mask PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${MASK_LIBRARY_DIR})