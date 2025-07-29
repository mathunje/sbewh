cmake_minimum_required(VERSION 3.28)

macro(SET_OPTION option value)
    set(${option} ${value} CACHE STRING INTERNAL FORCE)
endmacro()

project(SBE-Wannier-Houston)

set(SBE_WH_VERSION_MAJOR 0)
set(SBE_WH_VERSION_MINOR 1)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -O1" )
set(CMAKE_CXX_FLAGS "-O1")

# Global CUDA CXX flags/options
set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_EXTENSIONS OFF)
set(CMAKE_CUDA_FLAGS "${CMAKE_CXX_FLAGS} -g")

if ( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-attributes")
endif()


# add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")
option(SBE_WH_OPENMP "enable OpenMP" ON)
option(SBE_WH_MPI "enable MPI" ON)
option(SBE_WH_CUDA "enable Cuda" ON)

find_package(Boost CONFIG)


set(GIT_HASH_DIR gitHash)
include(${GIT_HASH_DIR}/CheckGit.cmake)
CheckGitSetup()

set(INCLUDE_DIR h/ ${GIT_HASH_DIR} util)
set(SRC_DIR src)
set(sbewhLibs git_version)


find_package(ZLIB REQUIRED)
list(APPEND sbewhLibs ${ZLIB_LIBRARIES})
list(APPEND INCLUDE_DIR ${ZLIB_INCLUDE_DIRS})


if ( SBE_WH_MPI )
    find_package(MPI REQUIRED)
    set_option(SBE_WH_MPI ON)
    list(APPEND INCLUDE_DIR ${MPI_INCLUDE_PATH})
    list(APPEND sbewhLibs MPI::MPI_CXX)
else()
    set_option(SBE_WH_MPI OFF)
endif()

if ( SBE_WH_OPENMP )
    find_package(OpenMP REQUIRED)
    set_option(SFP_OPENMP ON)
    list(APPEND sbewhLibs OpenMP::OpenMP_CXX)
else()
    set_option(SFP_OPENMP OFF)
endif()



find_package(LAPACK)
if ( LAPACK_FOUND )
    set_option(SBE_WH_LAPACK ON)
    set_option(SFP_LAPACK ON)
    list(APPEND sbewhLibs LAPACK::LAPACK)
else()
    set_option(SBE_WH_LAPACK OFF)
    set_option(SFP_LAPACK OFF)
endif()

if (SBE_WH_CUDA)
    find_package(CUDAToolkit COMPONENTS cublas cusolver cufft cudart)
    if( CMAKE_CUDA_ARCHITECTURES LESS 50)
        set(CMAKE_CUDA_ARCHITECTURES 50 60 70 75 80 86)
    endif()
    enable_language(CUDA)

    set(cudaCPP
        ${SRC_DIR}/CuDiagonalizator.cpp
        ${SRC_DIR}/CuWhFourierTransform.cpp
        ${SRC_DIR}/CuPropagator.cpp
        ${SRC_DIR}/CuObserver.cpp
        ${SRC_DIR}/runSBE.cpp
       )
    set_source_files_properties(${cudaCPP} PROPERTIES LANGUAGE CUDA)
    list(APPEND INCLUDE_DIR ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    list(APPEND sbewhLibs CUDA::cublas CUDA::cusolver CUDA::cufft CUDA::cudart)
else()

endif()


include_directories(${INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${PROJECT_BINARY_DIR} )


configure_file(sbewhConfig.h.in sbewhConfig.h)


add_executable(sbewh
               ${SRC_DIR}/main.cpp
               ${SRC_DIR}/runWHtransform.cpp
               ${SRC_DIR}/runFermiLevel.cpp
               ${SRC_DIR}/runSBE.cpp
               ${SRC_DIR}/runPulse.cpp
               #
               ${SRC_DIR}/Parameter.cpp
               ${SRC_DIR}/io/GeneralParameter.cpp
               ${SRC_DIR}/io/ParallelParameter.cpp
               ${SRC_DIR}/io/PulseParameter.cpp
               ${SRC_DIR}/io/TightBindingParameter.cpp
               ${SRC_DIR}/io/FourierTransformParameter.cpp
               ${SRC_DIR}/io/DiagonalizationParameter.cpp
               ${SRC_DIR}/io/PropagationParameter.cpp
               ${SRC_DIR}/io/KspaceRegionParameter.cpp
               ${SRC_DIR}/io/OutputParameter.cpp
               ${SRC_DIR}/io/util.cpp
               ${SRC_DIR}/io/save.cpp
               #
               ${SRC_DIR}/WhFourierTransform.cpp
               ${SRC_DIR}/StateContext.cpp
               ${SRC_DIR}/PropagatorBase.cpp
               ${SRC_DIR}/PropagatorCPU.cpp
               ${SRC_DIR}/ObserverBase.cpp
               ${SRC_DIR}/ObserverCPU.cpp
               #
               ${cudaCPP}
               #
               ${SRC_DIR}/mpiUtil.cpp

            external/cnpy/cnpy.cpp
            #
            util/src/Logger.cpp
            #
            #src/xyzIO.cpp
            #src/ObjWriter.cpp
            #
            #src/SparseRateOperator.cpp
            #src/PartialReductor.cpp
            #src/IndexedSparseMatrixVectorMultiplication.cpp
            #src/IndexedScalarProduct.cpp
            #src/RateExpectationValue.cpp
            #
            util/src/ParameterIO.cpp
            util/src/parameter/ParameterBase.cpp
            util/src/parameter/ParameterBool.cpp
            util/src/parameter/ParameterNumeric.cpp
            util/src/parameter/ParameterString.cpp
            util/src/parameter/ParameterMultiConfig.cpp
            util/src/detail/parse.cpp
            #
            util/src/pulse/Pulse1D.cpp
            util/src/pulse/GaussianPulse1D.cpp
            util/src/pulse/PureGauss1D.cpp
            util/src/pulse/Sin2Pulse1D.cpp
            util/src/pulse/Sin2RampPulse1D.cpp
            util/src/pulse/Pulse3D.cpp
            #
            util/src/w90util.cpp
            #
            util/src/DiagonalizatorSingle.cpp
            util/src/DiagonalizatorMultiple.cpp
            )


target_include_directories(sbewh PUBLIC h util/h util/h/pulse util/h/parameter util/h/detail external/cnpy )

target_link_libraries(sbewh ${sbewhLibs})
