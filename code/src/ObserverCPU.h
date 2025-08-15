#ifndef SBE_WH_OBSERVER_CPU_H
#define SBE_WH_OBSERVER_CPU_H

#include<algorithm>
#include<cmath>

#include "sbewhConfig.h"
#include "ObserverBase.h"
#include "PropagatorCPU.h"

#include "util/DiagonalizatorSingle.h"

#ifdef SBE_WH_MPI
#include <mpi.h>
#include "parallel/mpiUtil.h"
#endif

class ObserverCPU : public ObserverBase {
    WhFourierTransform &wft;
    HkFullVec fftMeshOps;
    long unsigned index;
    bool localMerging;

#ifdef SBE_WH_OPENMP
    std::vector< DiagonalizatorSingle > diag;
    std::vector< std::vector<std::complex<double>> > temp;
    std::vector< std::vector<double> > memEigenValues;
    std::vector< std::vector<double> > memOccWeights;
    std::vector< std::vector<double> > memRegionWeights;
#else
    DiagonalizatorSingle diag;
    std::vector<std::complex<double>> temp;
    std::vector<double> memEigenValues;
    std::vector<double> memOccWeights;
    std::vector<double> memRegionWeights;
#endif
    double calcWeight(unsigned kpt, const std::vector<GeomVector3d> &fineShifts, const ObserverRegion_t &region, const GeomVector3d &A);
    void mpi_reduceExpectationValuesEntry(ExpectationValuesEntry_t &ev);
public:
    ObserverCPU(const StateContext &sc, WhFourierTransform &wft);
    void reset();
    void setLocalMerging(bool m) { localMerging = m; }
    void operator()(const CPUstate &s, double t);
    void mpi_reduceExpectationValues();
    void normalizeExpectationValues();
};

#endif
