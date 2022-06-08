add_library(Mask STATIC)
target_include_directories(Mask 
    PUBLIC ${MASK_INCLUDE_DIR}
)

target_sources(Mask PRIVATE
    AngularDistribution.cpp
    DecaySystem.cpp
    EnergyLoss.cpp
    LayeredTarget.cpp
    LegendrePoly.cpp
    MaskApp.cpp
    MaskFile.cpp
    MassLookup.cpp
    Nucleus.cpp
    OneStepSystem.cpp
    RandomGenerator.cpp
    Reaction.cpp
    ReactionSystem.cpp
    Rotation.cpp
    Stopwatch.cpp
    Target.cpp
    ThreeStepSystem.cpp
    TwoStepSystem.cpp
    Vec3.cpp
    Vec4.cpp
)

set_target_properties(Mask PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${MASK_LIBRARY_DIR})