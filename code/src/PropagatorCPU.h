#ifndef SBE_WH_PROPAGATOR_CPU_H
#define SBE_WH_PROPAGATOR_CPU_H

#include "sbewhConfig.h"
#include "util/DiagonalizatorSingle.h"
#include "WhFourierTransform.h"
#include "PropagatorBase.h"

typedef std::vector<std::complex<double>> CPUstate;

#ifdef SBE_WH_LAPACK
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 double
#include<mkl_lapacke.h>
#include <mkl.h>
#endif

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
    inline void dephaseRhoH(std::complex<double> *rhoH, const double *energy);
    inline void dephaseRhoH_upper(std::complex<double> *rhoH, const double *energy);
    inline void manualCoherent(std::complex<double> *dsdtk, const std::complex<double> *Hk, const std::complex<double> *rhoK);
    inline void addManualDephasing(std::complex<double> *dsdtk, const std::complex<double> *rhoK,
                                   const std::complex<double> *U, const double* energy,
                                   std::complex<double> * rhoH, std::complex<double> *T);
#ifdef SBE_WH_LAPACK
    inline void lapackCoherent(std::complex<double> *dsdtk, const std::complex<double> *Hk, const std::complex<double> *rhoK);

    inline void addLapackDephasing(std::complex<double> *dsdtk, const std::complex<double> *rhoK,
                                   const std::complex<double> *U, const double* energy,
                                   std::complex<double> * rhoH, std::complex<double> *T);
#endif
public:
    PropagatorCPU(const StateContext &sc, WhFourierTransform &wft);
    void reset();
    CPUstate calcInitialState(double fermiLevel);
    void operator() (const CPUstate &state, CPUstate &dsdt, const double t);
};

#endif
