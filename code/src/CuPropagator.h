#ifndef SBE_WH_PROPAGATOR_GPU_H
#define SBE_WH_PROPAGATOR_GPU_H

#include "CuWhFourierTransform.h"
#include "CuDiagonalizator.h"
#include "PropagatorBase.h"

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include "cublas_v2.h"
#include "cuUtils.h"


// 2 consecutive doubles are interpretes as one complex number
typedef thrust::device_vector< double >  GPUstate;


class CuPropagator : public PropagatorBase {
    CuWhFourierTransformHk &wft;
    CuDiagonalizator &cud;
    // not owned
    cudaStream_t stream;
    cuDoubleComplex * d_H;
    // self managed
    cublasHandle_t cublasHandle;
    cuDoubleComplex * d_intermediate;
    cuDoubleComplex * d_rhoH;
    std::vector< std::complex<double> > specialValues_h;
    std::vector< std::complex<double> > specialValues_dsdt_h;
public:
    CuPropagator(cudaStream_t stream, const StateContext &sc, CuWhFourierTransformHk &wft, cuDoubleComplex * d_H, CuDiagonalizator &cud);
    ~CuPropagator();
    GPUstate calcInitialState(double fermiLevel);
    void operator() (const GPUstate &state, GPUstate &dsdt, const double t);
};

#endif
