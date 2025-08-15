target_sources(sbewh PRIVATE
    DiagonalizationParameter.cpp
    FourierTransformParameter.cpp
    GeneralParameter.cpp
    KspaceRegionParameter.cpp
    OutputParameter.cpp
    ParallelParameter.cpp
    ParameterIO.cpp
    PropagationParameter.cpp
    PulseParameter.cpp
    TightBindingParameter.cpp
    save.cpp
    util.cpp
)

add_subdirectory(parameter)
add_subdirectory(detail)
