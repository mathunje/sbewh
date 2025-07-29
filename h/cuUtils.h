#ifndef SBE_WH_CU_UTILS_H
#define SBE_WH_CU_UTILS_H


#include<string>

#include <cufft.h>
static std::string cufftGetErrorString(cufftResult error)
{
    switch (error)
    {
        case CUFFT_SUCCESS:
            return "CUFFT_SUCCESS";
        case CUFFT_INVALID_PLAN:
            return "CUFFT_INVALID_PLAN";
        case CUFFT_ALLOC_FAILED:
            return "CUFFT_ALLOC_FAILED";
        case CUFFT_INVALID_TYPE:
            return "CUFFT_INVALID_TYPE";
        case CUFFT_INVALID_VALUE:
            return "CUFFT_INVALID_VALUE";
        case CUFFT_INTERNAL_ERROR:
            return "CUFFT_INTERNAL_ERROR";
        case CUFFT_EXEC_FAILED:
            return "CUFFT_EXEC_FAILED";
        case CUFFT_SETUP_FAILED:
            return "CUFFT_SETUP_FAILED";
        case CUFFT_INVALID_SIZE:
            return "CUFFT_INVALID_SIZE";
        case CUFFT_UNALIGNED_DATA :
            return "CUFFT_UNALIGNED_DATA";
        case CUFFT_INCOMPLETE_PARAMETER_LIST:
            return "CUFFT_INCOMPLETE_PARAMETER_LIST";
        case CUFFT_INVALID_DEVICE:
            return "CUFFT_INVALID_DEVICE";
        case CUFFT_PARSE_ERROR:
            return "CUFFT_PARSE_ERROR";
        case CUFFT_NO_WORKSPACE:
            return "CUFFT_NO_WORKSPACE";
        case CUFFT_NOT_IMPLEMENTED:
            return "CUFFT_NOT_IMPLEMENTED";
        case CUFFT_LICENSE_ERROR:
            return "CUFFT_LICENSE_ERROR";
        case CUFFT_NOT_SUPPORTED:
            return "CUFFT_NOT_SUPPORTED";
    }
    return "unkown (" + std::to_string(static_cast<int>(error)) + ")";
}

#define CUDA_CHECK(err)                                                                            \
    do {                                                                                           \
        cudaError_t err_ = (err);                                                                  \
        if (err_ != cudaSuccess) {                                                                 \
            fprintf(stderr, "CUDA error %s at %s:%d\n",                                            \
                    cudaGetErrorString(err_), __FILE__, __LINE__);                                 \
        }                                                                                          \
    } while (0)

#define CUSOLVER_CHECK(err)                                                                        \
    do {                                                                                           \
        cusolverStatus_t err_ = (err);                                                             \
        if (err_ != CUSOLVER_STATUS_SUCCESS) {                                                     \
            fprintf(stderr, "cusolver error %d at %s:%d\n", err_, __FILE__, __LINE__);             \
        }                                                                                          \
    } while (0)




#define CUFFT_CHECK(err)                                                                           \
    do {                                                                                           \
        cufftResult err_ = static_cast<cufftResult>(err);                                          \
        if ( err_ != CUFFT_SUCCESS ) {                                                             \
            fprintf(stderr, "cufft error %s at %s:%d\n", cufftGetErrorString(err_).c_str()         \
                    , __FILE__, __LINE__);                                                         \
        }                                                                                          \
    } while (0)


#define CUBLAS_CHECK(err)                                                                          \
    do {                                                                                           \
        cublasStatus_t err_ = (err);                                                               \
        if (err_ != CUBLAS_STATUS_SUCCESS) {                                                       \
            fprintf(stderr, "cublas error %s at %s:%d\n", cublasGetStatusName(err_),               \
                    __FILE__, __LINE__);                                                           \
        }                                                                                          \
    } while (0)


#endif
