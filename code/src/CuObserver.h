#ifndef SBE_WH_CU_OBSERVER_H
#define SBE_WH_CU_OBSERVER_H

#include "ObserverBase.h"
#include "CuPropagator.h"


struct ExpectationIndexMapper{
private:
    unsigned cid = 0;
    unsigned nextId(unsigned cnt=1) { unsigned r = cid; cid += cnt; return r; }
public:
    const unsigned weights;
    const unsigned j;
    const unsigned occWan;
    const unsigned occHam;
    const unsigned numOps;

    ExpectationIndexMapper(const StateContext &sc) :
        weights(nextId() ),
        j(nextId(sc.dim) ),
        occWan(nextId(sc.numWann)),
        occHam(nextId(sc.numWann)),
        numOps(nextId(0))
        {}
};


class CuObserver : public ObserverBase {
    CuWhFourierTransformHkFull &wft;
    CuDiagonalizator &cud;
    cudaStream_t stream;
    cuDoubleComplex * d_H;
    // self managed
    ExpectationIndexMapper eim;
    cublasHandle_t cublasHandle;
    cuDoubleComplex * d_temp1;
    cuDoubleComplex * d_temp2;
    double * d_kBaseShifts;
    ObserverRegion_t * d_observerRegion;
    double * d_weights;
    cuDoubleComplex * d_complexWeights;
    cuDoubleComplex * d_wan;
    double * d_kexp;
    double * d_regionOps;
    std::vector<double> regionOps;
    std::vector<std::complex<double>> wan;
    std::vector< std::complex<double> > specialValues_h;
    unsigned numOps = 0;
public:
    CuObserver(cudaStream_t stream, const StateContext &sc, CuWhFourierTransformHkFull &wft, cuDoubleComplex * d_H, CuDiagonalizator &cud);
    ~CuObserver();
    void operator()(const GPUstate &s, double t);
};

#endif
