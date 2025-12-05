#include "PropagatorCPU.h"

PropagatorCPU::PropagatorCPU(const StateContext &sc, WhFourierTransform &wft) :
    PropagatorBase(sc), wft(wft), fftMeshOps( wft.createHkVec() ),
    tempStride( ((sc.numWann * sc.numWann + alignmentFct -1) / alignmentFct) * alignmentFct)
#ifndef SBE_WH_OPENMP
    , diag(sc.param.diag.mode, sc.numWann), memEigenValues(sc.numWann), temp(2 * tempStride)
#endif
{
#ifdef SBE_WH_OPENMP
    memEigenValues.resize(omp_get_max_threads());
    temp.resize(omp_get_max_threads());
    for(unsigned tid=0; tid<omp_get_max_threads(); tid++){
        DiagonalizatorSingle d(sc.param.diag.mode, sc.numWann);
        d.setSweepCount(sc.param.diag.maxSweepCount);
        diag.push_back(d);
        memEigenValues[tid].resize(sc.numWann);
        temp[tid].resize(2*tempStride);
    }
#else
    diag.setSweepCount(sc.param.diag.maxSweepCount);
#endif
    reset();
}

void PropagatorCPU::reset()
{
    PropagatorBase::reset();
}

CPUstate PropagatorCPU::calcInitialState(double fermiLevel)
{
    GeomVector3d Ezero({0, 0, 0});
    const auto &shifts = wft.getFineMeshFracShifts(sc.kGlobalShift_au, sc.NrProcStart, sc.NrProcCount);
    wft.transformMeshFFT(fftMeshOps, shifts, Ezero);
    unsigned nKpts = fftMeshOps.getNumKpts();
    unsigned numWann = sc.numWann;
    CPUstate res(sc.si.specialIndexCount + sc.si.paddedMatSize * nKpts);
    for(unsigned dir=0; dir<3; dir++)
        res[sc.si.A[dir]] = 0;
#ifdef SBE_WH_OPENMP
    const unsigned occAlign = 32;
    std::vector<double> occupations( occAlign *((numWann+occAlign-1) / occAlign) * omp_get_max_threads(), 0);
#else
    std::vector<double> occupations(numWann, 0);
#endif
    #pragma omp parallel for schedule(static)
    for(unsigned kpt=0; kpt<fftMeshOps.getNumKpts(); kpt++){
        std::complex<double> * rhoK = res.data() + sc.si.rhoK(kpt);
#ifdef SBE_WH_OPENMP
        DiagonalizatorSingle &d = diag[omp_get_thread_num()];
        double * e = memEigenValues[omp_get_thread_num()].data();
#else
        DiagonalizatorSingle &d = diag;
        double * e = memEigenValues.data();
#endif
        std::complex<double> * U = fftMeshOps.H0(kpt);
        d.diag(U, e);
        for(unsigned i=0; i<numWann; i++)
            e[i] = 1 / ( 1 + exp( (e[i]-fermiLevel) / sc.param.tb.gsTemp) );
        for(unsigned a=0; a<numWann; a++)
            for(unsigned b=0; b<numWann; b++){
                std::complex<double> sum = 0;
                for(unsigned i=0; i<numWann; i++)
                    sum += U[i*numWann+b] * e[i] * std::conj(U[i*numWann+a]);
                rhoK[a*numWann+b] = sum;
            }
    }
    return res;
}


inline void PropagatorCPU::dephaseRhoH(std::complex<double> *rhoH, const double *energy)
{
    unsigned numWann = sc.numWann;
    if ( sc.param.prop.soothingWidth < 0 ){
        for(unsigned a=0; a<numWann; a++)
            for(unsigned b=0; b<numWann; b++)
                rhoH[a*numWann+b] = - (a!=b) * sc.param.prop.relaxationF2 * rhoH[a*numWann+b];
    } else {
        for(unsigned a=0; a<numWann; a++)
            for(unsigned b=0; b<numWann; b++){
                double dE = (energy[a] - energy[b]) / sc.param.prop.soothingWidth;
                rhoH[a*numWann+b] = -(1 - std::exp(-dE*dE) ) * sc.param.prop.relaxationF2 * rhoH[a*numWann+b];
            }
    }
}
inline void PropagatorCPU::dephaseRhoH_upper(std::complex<double> *rhoH, const double *energy)
{
    unsigned numWann = sc.numWann;
    if ( sc.param.prop.soothingWidth < 0 ){
        for(unsigned a=0; a<numWann; a++){
            rhoH[a*numWann+a] = 0;
            for(unsigned b=a+1; b<numWann; b++)
                rhoH[a*numWann+b] = -sc.param.prop.relaxationF2 * rhoH[a*numWann+b];
        }
    } else {
        for(unsigned a=0; a<numWann; a++){
            rhoH[a*numWann+a] = 0;
            for(unsigned b=a+1; b<numWann; b++){
                double dE = (energy[a] - energy[b]) / sc.param.prop.soothingWidth;
                rhoH[a*numWann+b] = -(1 - std::exp(-dE*dE) ) * sc.param.prop.relaxationF2 * rhoH[a*numWann+b];
            }
        }
    }
}


inline void PropagatorCPU::manualCoherent(std::complex<double> *dsdtk, const std::complex<double> *Hk, const std::complex<double> *rhoK)
{
    unsigned numWann = sc.numWann;
    for(unsigned a=0; a<numWann; a++)
        for(unsigned b=0; b<numWann; b++){
            std::complex<double> sum = 0;
            for(unsigned i=0; i<numWann; i++)
                sum += Hk[a*numWann+i] * rhoK[i*numWann+b] - Hk[i*numWann+b] * rhoK[a*numWann+i];
            dsdtk[a*numWann+b] = std::complex<double>(0, -1) * sum;
        }
}

inline void PropagatorCPU::addManualDephasing(std::complex<double> *dsdtk, const std::complex<double> *rhoK,
                                              const std::complex<double> *U, const double* energy,
                                              std::complex<double> * rhoH, std::complex<double> *T)
{
    unsigned numWann = sc.numWann;
    for(unsigned a=0; a<numWann; a++)
        for(unsigned b=0; b<numWann; b++){
            std::complex<double> sum = 0;
            for(unsigned i=0; i<numWann; i++)
                sum += U[a*numWann+i] * rhoK[i*numWann+b];
            T[a*numWann+b] = sum;
        }
    for(unsigned a=0; a<numWann; a++)
        for(unsigned b=0; b<numWann; b++){
            std::complex<double> sum = 0;
            for(unsigned i=0; i<numWann; i++)
                sum += T[a*numWann+i] * std::conj(U[b*numWann+i]);
            rhoH[a*numWann+b] = sum;
        }
    dephaseRhoH(rhoH, energy);
    for(unsigned a=0; a<numWann; a++)
        for(unsigned b=0; b<numWann; b++){
            std::complex<double> sum = 0;
            for(unsigned i=0; i<numWann; i++)
                sum += std::conj(U[i*numWann+a]) * rhoH[i*numWann+b];
            T[a*numWann+b] = sum;
        }
    for(unsigned a=0; a<numWann; a++)
        for(unsigned b=0; b<numWann; b++){
            std::complex<double> sum = 0;
            for(unsigned i=0; i<numWann; i++)
                sum += T[a*numWann+i] * U[i*numWann+b];
            dsdtk[a*numWann+b] += sum;
        }
}

#ifdef SBE_WH_LAPACK
inline void PropagatorCPU::lapackCoherent(std::complex<double> *dsdtk, const std::complex<double> *Hk, const std::complex<double> *rhoK)
{
    int numWann = (int)sc.numWann;
    char S = 'L';
    char U = 'L';
    std::complex<double> alpha(0, 1);
    std::complex<double> beta(0, 0);
    zhemm(&S, &U, &numWann, &numWann, &alpha, Hk, &numWann,
          rhoK, &numWann, &beta,
          dsdtk, &numWann);
    S = 'R';
    alpha = std::complex<double>(0, -1);
    beta = std::complex<double>(1, 0);
    zhemm(&S, &U, &numWann, &numWann, &alpha, Hk, &numWann,
          rhoK, &numWann, &beta,
          dsdtk, &numWann);
}

inline void PropagatorCPU::addLapackDephasing(std::complex<double> *dsdtk, const std::complex<double> *rhoK,
                                              const std::complex<double> *U, const double* energy,
                                              std::complex<double> * rhoH, std::complex<double> *T)
{
    int numWann = (int)sc.numWann;
    char S = 'L';
    char Uplo = 'L';
    std::complex<double> alpha(1, 0);
    std::complex<double> beta(0, 0);
    zhemm(&S, &Uplo, &numWann, &numWann, &alpha, rhoK, &numWann,
          U, &numWann, &beta,
          T, &numWann);
    char TA = 'C', TB = 'N';
    zgemmt(&Uplo, &TA, &TB, &numWann, &numWann, &alpha, T, &numWann, U, &numWann, &beta, rhoH, &numWann);
    dephaseRhoH_upper(rhoH, energy);
    S = 'R';
    zhemm(&S, &Uplo, &numWann, &numWann, &alpha, rhoH, &numWann,
          U, &numWann, &beta,
          T, &numWann);
    std::complex<double> * dephTmp = rhoH;
    zgemmt(&Uplo, &TB, &TA, &numWann, &numWann, &alpha, T, &numWann, U, &numWann, &beta, dephTmp, &numWann);
    for(int a=0; a<numWann; a++){
        dsdtk[a*(numWann+1)] += std::real(dephTmp[a*(numWann+1)]);
        for(int b=a+1; b<numWann; b++){
            std::complex<double> tmp = dephTmp[a*numWann+b];
            dsdtk[a*numWann+b] += tmp;
            dsdtk[b*numWann+a] += std::conj(tmp);
        }
    }
}
#endif

void PropagatorCPU::operator() (const CPUstate &state, CPUstate &dsdt, const double t)
{
    auto startTime = std::chrono::steady_clock::now();
    std::fill(dsdt.begin(), dsdt.end(), 0);
    unsigned numWann = sc.numWann;
    GeomVector3d E = sc.getE(t);
    GeomVector3d A({ 0, 0, 0 });
    for(unsigned i=0; i<sc.dim; i++){
        A[i] = std::real(state[sc.si.A[i]]);
        dsdt[sc.si.A[i]] = -E[i];
    }
    GeomVector3d kOff = sc.kGlobalShift_au + A;
    const auto &shifts = wft.getFineMeshFracShifts(kOff, sc.NrProcStart, sc.NrProcCount);
    wft.transformMeshFFT(fftMeshOps, shifts, E);
    #pragma omp parallel for schedule(static)
    for(unsigned kpt=0; kpt<fftMeshOps.getNumKpts(); kpt++){
        const std::complex<double> * rhoK = state.data() + sc.si.rhoK(kpt);
        const std::complex<double> * Hk = fftMeshOps.HE(kpt);
        std::complex<double> * dsdtk = dsdt.data() + sc.si.rhoK(kpt);
#ifdef SBE_WH_LAPACK
        lapackCoherent(dsdtk, Hk, rhoK);
#else
        manualCoherent(dsdt.data() + sc.si.rhoK(kpt), Hk, rhoK);
#endif
        // dephasing in Hamiltonian gauge
        if ( sc.param.prop.relaxationF2 > 0 ){
#ifdef SBE_WH_OPENMP
            DiagonalizatorSingle &d = diag[omp_get_thread_num()];
            double * energy = memEigenValues[omp_get_thread_num()].data();
            std::complex<double> * T = temp[omp_get_thread_num()].data();
            std::complex<double> * rhoH = temp[omp_get_thread_num()].data() + tempStride;
#else
            DiagonalizatorSingle &d = diag;
            double * energy = memEigenValues.data();
            std::complex<double> * T = temp.data();
            std::complex<double> * rhoH = temp.data() + tempStride;
#endif
            std::complex<double> * U = fftMeshOps.H0(kpt);
            d.diag(U, energy);
#ifdef SBE_WH_LAPACK
            addLapackDephasing(dsdtk, rhoK, U, energy, rhoH, T);
#else
            addManualDephasing(dsdtk, rhoK, U, energy, rhoH, T);
#endif
        }
    }
    if ( sc.param.general.timeIt )
        timeAndPrint(startTime);
}
