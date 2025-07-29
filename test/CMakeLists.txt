cmake_minimum_required(VERSION 3.18)

macro(SET_OPTION option value)
    set(${option} ${value} CACHE STRING INTERNAL FORCE)
endmacro()

project(SBE-Test)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -O1" )

if ( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-attributes")
endif()

set(INCLUDE_DIR ../h ../util/h ../util/h/parameter ../util/h/pulse ../util/h/detail)
set(SRC_DIR ../src)
set(TEST_DIR .)

find_package(GTest CONFIG REQUIRED)
set(testLinkLibs gtest gtest_main)

find_package(OpenMP)
if ( OpenMP_FOUND )
    set_option(SBE_WH_OPENMP ON)
    list(APPEND testLinkLibs OpenMP::OpenMP_CXX)
else()
    set_option(SBE_WH_OPENMP OFF)
endif()

find_package(LAPACK)
if ( LAPACK_FOUND )
    set_option(SBE_WH_LAPACK ON)
    list(APPEND testLinkLibs LAPACK::LAPACK)
else()
    set_option(SBE_WH_LAPACK OFF)
endif()

find_package(CUDAToolkit REQUIRED)
if (CUDAToolkit_FOUND)
    set_option(SBE_WH_CUDA ON)
    enable_language(CUDA)
    if( CMAKE_CUDA_ARCHITECTURES LESS 60)
        set(CMAKE_CUDA_ARCHITECTURES 60 70 75 80 86)
    endif()
    # set_source_files_properties(${INCLUDE_DIR}/CuDiagonalizator.h PROPERTIES LANGUAGE CUDA)
    list(APPEND INCLUDE_DIR ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    list(APPEND testLinkLibs CUDA::cusolver)
    set(cudaCPP ${TEST_DIR}/cuDiagonalizator.cpp ${SRC_DIR}/CuDiagonalizator.cpp)
    set_source_files_properties(${cudaCPP} PROPERTIES LANGUAGE CUDA)
else()
    set_option(SBE_WH_CUDA OFF)
endif()

configure_file(../sbewhConfig.h.in sbewhConfig.h)

include_directories(${INCLUDE_DIR} ${PROJECT_BINARY_DIR})

set(PULSE_DIR ../util/src/pulse)
set(PULSE_SRC ${PULSE_DIR}/Pulse1D.cpp
             ${PULSE_DIR}/PureGauss1D.cpp
             ${PULSE_DIR}/GaussianPulse1D.cpp
             ${PULSE_DIR}/Sin2Pulse1D.cpp
             ${PULSE_DIR}/Sin2RampPulse1D.cpp
             ${PULSE_DIR}/Pulse3D.cpp)

set(DETAIL_DIR ../util/src/detail)

set(PARAM_DIR ../util/src/parameter)
set(PARAM_SRC ${PARAM_DIR}/ParameterBase.cpp
              ${PARAM_DIR}/ParameterBool.cpp
              ${PARAM_DIR}/ParameterNumeric.cpp
              ${PARAM_DIR}/ParameterString.cpp
              ${PARAM_DIR}/ParameterMultiConfig.cpp)


add_executable(t ${TEST_DIR}/unitConversion.cpp
                 ${TEST_DIR}/geomVector.cpp
                 ${TEST_DIR}/parameter.cpp ${PARAM_SRC}
                 ${TEST_DIR}/parameterIO.cpp ../util/src/ParameterIO.cpp ${DETAIL_DIR}/parse.cpp
                 ${TEST_DIR}/pulse.cpp ${PULSE_SRC}
                 ${TEST_DIR}/hermiteSpline.cpp
                 ${TEST_DIR}/w90util.cpp ../util/src/w90util.cpp
                 ${cudaCPP})
target_link_libraries(t ${testLinkLibs})
