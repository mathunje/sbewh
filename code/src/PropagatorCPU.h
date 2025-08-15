#ifndef SBE_WH_PROPAGATOR_CPU_H
#define SBE_WH_PROPAGATOR_CPU_H

#include "sbewhConfig.h"
#include "util/DiagonalizatorSingle.h"
#include "WhFourierTransform.h"
#include "PropagatorBase.h"

typedef std::vector<std::complex<double>> CPUstate;

#ifdef SBE_WH_OPENMP
#include <omp.h>
#endif


class PropagatorCPU : public PropagatorBase {
    WhFourierTransform &wft;
    HkVec fftMeshOps;
    const unsigned alignmentFct = 64;
    const unsigned tempStride;
#ifdef SBE_WH_OPENMP
    std::vector<DiagonalizatorSingle> diag;
    std::vector< std::vector<double> > memEigenValues;
    std::vector< std::vector<std::complex<double>> > temp;
#else
    DiagonalizatorSingle diag;
    std::vector<double> memEigenValues;
    std::vector<std::complex<double>> temp;
#endif
    CPUstate calcInitialStateFromOps(double fermiLevel);
public:
    PropagatorCPU(const StateContext &sc, WhFourierTransform &wft);
    void reset();
    CPUstate calcInitialState(double fermiLevel);
    void operator() (const CPUstate &state, CPUstate &dsdt, const double t);
};

#endif
