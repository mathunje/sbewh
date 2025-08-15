cmake_minimum_required(VERSION 3.28)

macro(SET_OPTION option value)
    set(${option} ${value} CACHE STRING INTERNAL FORCE)
endmacro()

project(SBE-Wannier-Houston)

set(SBE_WH_VERSION_MAJOR 0)
set(SBE_WH_VERSION_MINOR 2)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -O1" )
set(CMAKE_CXX_FLAGS "-O1")

if ( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-attributes")
endif()

option(SBE_WH_BUILD_TESTS "compilation of tests" ON)

option(SBE_WH_OPENMP "enable OpenMP" OFF)
option(SBE_WH_MPI "enable MPI" OFF)
option(SBE_WH_CUDA "enable Cuda" OFF)
option(SBE_WH_LAPACK "enable LAPACK" OFF)


find_package(Boost CONFIG 1.53 REQUIRED)

set(sbewhLibs fftw3)
set(sbewhIncludeDirs "${PROJECT_BINARY_DIR}" "${PROJECT_BINARY_DIR}/generated")

if ( SBE_WH_OPENMP )
    find_package(OpenMP REQUIRED)
    list(APPEND sbewhLibs OpenMP::OpenMP_CXX)
endif()

if ( SBE_WH_MPI )
    find_package(MPI REQUIRED)
    list(APPEND sbewhIncludeDirs ${MPI_INCLUDE_PATH})
    list(APPEND sbewhLibs MPI::MPI_CXX)
endif()

if ( SBE_WH_LAPACK )
    find_package(LAPACK)
    list(APPEND sbewhLibs LAPACK::LAPACK)
endif()

if ( SBE_WH_CUDA )
    find_package(CUDAToolkit COMPONENTS cublas cusolver cufft cudart)
    enable_language(CUDA)
    set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_EXTENSIONS OFF)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CXX_FLAGS} -g")

    list(APPEND sbewhIncludeDirs ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    list(APPEND sbewhLibs CUDA::cublas CUDA::cusolver CUDA::cufft CUDA::cudart)
endif()

configure_file(sbewhConfig.h.in sbewhConfig.h)

include(cmake_modules/CheckGit.cmake)
CheckGitSetup()

add_executable(sbewh)
target_link_libraries(sbewh PRIVATE "${sbewhLibs}")
target_include_directories(sbewh PRIVATE "${sbewhIncludeDirs}")
target_sources(sbewh PRIVATE ${PROJECT_BINARY_DIR}/generated/git_version.cpp)


add_subdirectory(src)
if ( SBE_WH_BUILD_TESTS )
    add_executable(test_sbewh)
    target_link_libraries(test_sbewh PRIVATE "${sbewhLibs}")
    target_include_directories(test_sbewh PRIVATE "${sbewhIncludeDirs}")
    add_subdirectory(tests)
endif()


include(GNUInstallDirs)
find_package(bash-completion)
if(BASH_COMPLETION_FOUND)
  message(STATUS
    "Using bash completion dir ${BASH_COMPLETION_COMPLETIONSDIR}")
else()
  set (BASH_COMPLETION_COMPLETIONSDIR "${CMAKE_INSTALL_DATADIR}/bash-completion/completions")
  message (STATUS
    "Using fallback bash completion dir ${BASH_COMPLETION_COMPLETIONSDIR}")
endif()

install(FILES bash-completion/sbewh DESTINATION
  ${BASH_COMPLETION_COMPLETIONSDIR})

install(TARGETS sbewh)
