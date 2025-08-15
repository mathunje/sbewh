target_sources(sbewh PRIVATE
   main.cpp
   runFermiLevel.cpp
   runPulse.cpp
   runSBE.cpp
   runWHtransform.cpp
)

if (SBE_WH_CUDA)
    target_sources(sbewh PRIVATE
        CuRunSBE.cu
    )
endif()
