set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -O1" )


find_package(GTest CONFIG REQUIRED)
target_link_libraries(test_sbewh PRIVATE gtest gtest_main)



if (SBE_WH_CUDA)
    target_sources(test_sbewh PRIVATE
        cuDiagonalizator.cu ../src/CuDiagonalizator.cu
       )
endif()

set(PARAM_DIR ../src/io/parameter)
set(PARAM_SRC ${PARAM_DIR}/ParameterBase.cpp
              ${PARAM_DIR}/ParameterBool.cpp
              ${PARAM_DIR}/ParameterNumeric.cpp
              ${PARAM_DIR}/ParameterString.cpp
              ${PARAM_DIR}/ParameterMultiConfig.cpp)

set(PULSE_DIR ../src/pulse)
set(PULSE_SRC ${PULSE_DIR}/Pulse1D.cpp
             ${PULSE_DIR}/PureGauss1D.cpp
             ${PULSE_DIR}/GaussianPulse1D.cpp
             ${PULSE_DIR}/Sin2Pulse1D.cpp
             ${PULSE_DIR}/Sin2RampPulse1D.cpp
             ${PULSE_DIR}/Pulse3D.cpp)

target_include_directories(test_sbewh PRIVATE ${CMAKE_SOURCE_DIR}/src)

target_sources(test_sbewh PRIVATE
    unitConversion.cpp
    geomVector.cpp
    parameter.cpp ${PARAM_SRC}
    parameterIO.cpp ../src/io/ParameterIO.cpp ../src/io/detail/parse.cpp
    pulse.cpp ${PULSE_SRC}
    hermiteSpline.cpp
    w90util.cpp  ../src/util/w90util.cpp
)
