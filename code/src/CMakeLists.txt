target_sources(sbewh PRIVATE
    ObserverBase.cpp
    ObserverCPU.cpp
    Parameter.cpp
    PropagatorBase.cpp
    PropagatorCPU.cpp
    StateContext.cpp
    WhFourierTransform.cpp
)

if (SBE_WH_CUDA)
    target_sources(sbewh PRIVATE
        CuDiagonalizator.cu
        CuWhFourierTransform.cu
        CuPropagator.cu
        CuObserver.cu
       )
endif()

target_include_directories(sbewh PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(external)
add_subdirectory(io)
add_subdirectory(main)
add_subdirectory(parallel)
add_subdirectory(pulse)
add_subdirectory(util)
