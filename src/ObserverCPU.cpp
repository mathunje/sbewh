#include "ObserverCPU.h"

ExpectationValuesEntry_t ompExpectationValuesEntryInit(const ExpectationValuesEntry_t &rhs)
{
    ExpectationValuesEntry_t lhs;
    lhs.region.resize(rhs.region.size());
    lhs.occupationHamMax.resize(rhs.occupationHamMax.size(), std::numeric_limits<double>::lowest());
    lhs.occupationHamMin.resize(rhs.occupationHamMin.size(), std::numeric_limits<double>::max());
    for(unsigned rid=0; rid<rhs.region.size(); rid++){
        ExpectationValuesRegion_t &l = lhs.region[rid];
        const ExpectationValuesRegion_t &r = rhs.region[rid];
        l.kPointWeightSum = 0;
        for(unsigned dir=0; dir<3; dir++)
            l.j[dir] = 0;
        l.occupationHam.resize(r.occupationHam.size(), 0);
        l.occupationWan.resize(r.occupationWan.size(), 0);
        l.wan.resize(r.wan.size(), 0.0);
    }
    return lhs;
}

void ompExpectationValuesEntryReduce(ExpectationValuesEntry_t &lhs, const ExpectationValuesEntry_t &rhs)
{
    for(unsigned rid=0; rid<rhs.region.size(); rid++){
        ExpectationValuesRegion_t &l = lhs.region[rid];
        const ExpectationValuesRegion_t &r = rhs.region[rid];
        l.kPointWeightSum += r.kPointWeightSum;
        for(unsigned i=0; i<3; i++)
            l.j[i] += r.j[i];
        for(unsigned i=0; i<r.occupationHam.size(); i++){
            l.occupationHam[i] += r.occupationHam[i];
            l.occupationWan[i] += r.occupationWan[i];
        }
        for(unsigned i=0; i<r.wan.size(); i++)
            l.wan[i] += r.wan[i];
    }
    for(unsigned i=0; i<rhs.occupationHamMin.size(); i++){
        lhs.occupationHamMin[i] = std::min(lhs.occupationHamMin[i], rhs.occupationHamMin[i]);
        lhs.occupationHamMax[i] = std::max(lhs.occupationHamMax[i], rhs.occupationHamMax[i]);
    }
}


#pragma omp declare reduction( all : ExpectationValuesEntry_t : \
         ompExpectationValuesEntryReduce(omp_out, omp_in) ) \
         initializer(omp_priv = ompExpectationValuesEntryInit(omp_orig))


/* to reduce the number of MPI calls to one for each reduce type */
inline void appendPackedExpectationValuesEntry(const ExpectationValuesEntry_t &ev, std::vector<double> &sumEntries,
                                                    std::vector<double> &minEntries, std::vector<double> &maxEntries)
{
    for(unsigned rid=0; rid<ev.region.size(); rid++){
        const ExpectationValuesRegion_t &r = ev.region[rid];
        sumEntries.push_back(r.kPointWeightSum);
        for(unsigned i=0; i<3; i++)
            sumEntries.push_back(r.j[i]);
        for(unsigned i=0; i<r.occupationHam.size(); i++){
            sumEntries.push_back(r.occupationHam[i]);
            sumEntries.push_back(r.occupationWan[i]);
        }
        for(unsigned i=0; i<r.wan.size(); i++){
            sumEntries.push_back(std::real(r.wan[i]));
            sumEntries.push_back(std::imag(r.wan[i]));
        }
    }
    for(unsigned i=0; i<ev.occupationHamMin.size(); i++){
        minEntries.push_back(ev.occupationHamMin[i]);
        maxEntries.push_back(ev.occupationHamMax[i]);
    }
}


/* the {sum|min|max}I are starting indices to parse from and are advanced during parsing  */
inline void unpackExpectationValuesEntry(ExpectationValuesEntry_t &ev, const std::vector<double> &sumEntries, unsigned &sumI,
                                        const std::vector<double> & minEntries, unsigned &minI,
                                        const std::vector<double> & maxEntries, unsigned &maxI)
{
    for(unsigned rid=0; rid<ev.region.size(); rid++){
        ev.region[rid].kPointWeightSum = sumEntries[sumI++];
        for(unsigned i=0; i<3; i++)
            ev.region[rid].j[i] = sumEntries[sumI++];
        for(unsigned i=0; i<ev.region[rid].occupationHam.size(); i++){
            ev.region[rid].occupationHam[i] = sumEntries[sumI++];
            ev.region[rid].occupationWan[i] = sumEntries[sumI++];
        }
        for(unsigned i=0; i<ev.region[rid].wan.size(); i++){
            double re = sumEntries[sumI++];
            double im = sumEntries[sumI++];
            ev.region[rid].wan[i] = std::complex<double>(re, im);
        }
    }
    for(unsigned i=0; i<ev.occupationHamMin.size(); i++){
        ev.occupationHamMin[i] = minEntries[minI++];
        ev.occupationHamMax[i] = maxEntries[maxI++];
    }
}



/****************************************************************************************************/

ObserverCPU::ObserverCPU(const StateContext &sc, WhFourierTransform &wft): ObserverBase(sc),
    wft(wft), fftMeshOps( wft.createHkFullVec() )
#ifndef SBE_WH_OPENMP
    , diag(sc.param.diag.mode, sc.numWann), temp(sc.numWann * sc.numWann), memEigenValues(sc.numWann),
      memOccWeights(sc.numWann * sc.numWann), memRegionWeights(region.size())
#endif
{
#ifdef SBE_WH_OPENMP
    temp.resize(omp_get_max_threads());
    memEigenValues.resize(omp_get_max_threads());
    memOccWeights.resize(omp_get_max_threads());
    memRegionWeights.resize(omp_get_max_threads());
    for(unsigned tid=0; tid<omp_get_max_threads(); tid++){
        DiagonalizatorSingle d(sc.param.diag.mode, sc.numWann);
        d.setSweepCount(sc.param.diag.maxSweepCount);
        diag.push_back(d);
        memEigenValues[tid].resize(sc.numWann);
        memOccWeights[tid].resize(sc.numWann * sc.numWann);
        memRegionWeights[tid].resize(region.size());
        temp[tid].resize(sc.numWann * sc.numWann);
    }
#else
    diag.setSweepCount(sc.param.diag.maxSweepCount);
#endif
    reset();
}

void ObserverCPU::reset()
{
    index = 0;
    localMerging = false;
    ObserverBase::reset();
}

double ObserverCPU::calcWeight(unsigned kpt, const std::vector<GeomVector3d> &fineShifts, const ObserverRegion_t &region, const GeomVector3d &A)
{
    GeomVector3d pos = wft.getFFTfracShift(kpt, fineShifts, (region.movingFrame ? -A : GeomVector3d({0, 0, 0}) ) );
    pos -= region.k;
    double w = 1;
    for(unsigned dir=0; dir<3; dir++){
        double cPos = pos[dir] - std::round(pos[dir]);
        w *= std::exp( -cPos * cPos * region.invSigma2);
    }
    return w;
}

void ObserverCPU::operator()(const CPUstate &s, double t)
{
    basicLog(t);
    unsigned numWann = sc.numWann;
    unsigned maxDim = sc.dim;
    unsigned numKpts = fftMeshOps.getNumKpts();
    unsigned realIndex = localMerging ? index % expValues.size() : index;
    bool storeRegionDensity = getStoreRegionDensity( realIndex );
    if ( ! localMerging )
        expValues.push_back( prepareExpValue(t, storeRegionDensity) );
    ExpectationValuesEntry_t &ev = expValues[realIndex];
    GeomVector3d A;
    for(unsigned i=0; i<3; i++){
        A[i] = std::real(s[sc.si.A[i]]);
        ev.A[i] = A[i];
    }
    GeomVector3d kOff = sc.kGlobalShift_au + A;
    const auto &fineShifts = wft.getFineMeshFracShifts(kOff, sc.NrProcStart, sc.NrProcCount);
    wft.transformMeshFFT(fftMeshOps, fineShifts);
    # pragma omp parallel for schedule(static) reduction( all : ev )
    for(unsigned kpt=0; kpt<numKpts; kpt++){
#ifdef SBE_WH_OPENMP
        std::complex<double> * T = temp[omp_get_thread_num()].data();
        double * energy = memEigenValues[omp_get_thread_num()].data();
        double * occWeights = memOccWeights[omp_get_thread_num()].data();
        double * regionWeights = memRegionWeights[omp_get_thread_num()].data();
        DiagonalizatorSingle &d = diag[omp_get_thread_num()];
#else
        std::complex<double> * T = temp.data();
        double * energy = memEigenValues.data();
        double * regionWeights = memRegionWeights.data();
        DiagonalizatorSingle &d = diag;
#endif
        for(unsigned rid=0; rid<region.size(); rid++){
            double w = calcWeight(kpt, fineShifts, region[rid], A);
            regionWeights[rid] = w;
            ev.region[rid].kPointWeightSum += w;
        }
        const std::complex<double> * rhoK = s.data() + sc.si.rhoK(kpt);
        const std::complex<double> * Hk = fftMeshOps.H0(kpt);
        for(unsigned a=0; a<numWann; a++)
            for(unsigned b=0; b<numWann; b++){
                std::complex<double> sum = 0;
                for(unsigned i=0; i<numWann; i++)
                    sum += Hk[a*numWann+i] * rhoK[i*numWann+b] - Hk[i*numWann+b] * rhoK[a*numWann+i];
                T[a*numWann+b] = sum;
            }
        std::complex<double> jk[3] = {0, 0, 0};
        for(unsigned dir=0; dir<maxDim; dir++){
            const std::complex<double> * Dk = fftMeshOps.D(kpt, dir);
            const std::complex<double> * dH0dk = fftMeshOps.dH0dk(kpt, dir);
            for(unsigned a=0; a<numWann; a++)
                for(unsigned b=0; b<numWann; b++)
                    jk[dir] += std::complex<double>{0, 1} * Dk[a*numWann+b] * T[b*numWann+a]
                               - dH0dk[a*numWann+b] * rhoK[b*numWann+a];
        }
        for(unsigned rid=0; rid<region.size(); rid++){
            double weight = regionWeights[rid];
            for(unsigned dir=0; dir<3; dir++)
                ev.region[rid].j[dir] += std::real(jk[dir]) * weight;
        }
        std::complex<double> * U = fftMeshOps.H0(kpt);
        d.diag(U, energy);
        if ( sc.param.prop.occupationSmearingWidth >= 0 ){
            for(unsigned a=0; a<numWann; a++){
                double wSum = 0;
                for(unsigned b=0; b<numWann; b++){
                    double dE = (energy[a]-energy[b]) / sc.param.prop.occupationSmearingWidth;
                    double w = std::exp(-dE*dE);
                    occWeights[a*numWann+b] = w;
                    wSum += w;
                }
                for(unsigned b=0; b<numWann; b++)
                    occWeights[a*numWann+b] /= wSum;
            }
        }
        for(unsigned bi=0; bi<numWann; bi++){
            std::complex<double> sum = 0;
            for(unsigned a=0; a<numWann; a++)
                for(unsigned b=0; b<numWann; b++)
                    sum += std::conj(U[bi*numWann+b]) * rhoK[a*numWann+b] * U[bi*numWann+a];
            T[bi] = sum; // raw occupations
        }
        for(unsigned bi=0; bi<numWann; bi++){
            std::complex<double> weightedOcc = 0;
            if ( sc.param.prop.occupationSmearingWidth >= 0 ){
                for(unsigned i=0; i<numWann; i++)
                    weightedOcc += T[i] * occWeights[numWann*bi + i];
            } else {
                weightedOcc = T[bi];
            }
            ev.occupationHamMin[bi] = std::min(ev.occupationHamMin[bi], std::real(weightedOcc) );
            ev.occupationHamMax[bi] = std::max(ev.occupationHamMax[bi], std::real(weightedOcc) );
            for(unsigned rid=0; rid<region.size(); rid++){
                double weight = regionWeights[rid];
                ev.region[rid].occupationHam[bi] += std::real(weightedOcc) * weight;
                ev.region[rid].occupationWan[bi] += std::real(rhoK[bi*(numWann+1)]) * weight;
            }
        }
        if( storeRegionDensity ){
            for(unsigned rid=0; rid<region.size(); rid++){
                double weight = regionWeights[rid];
                for(unsigned i=0; i<numWann*numWann; i++)
                    ev.region[rid].wan[i] += rhoK[i] * weight;
            }
        }
    }
    if ( ! sc.param.par.independentMpiIntegrations )
        mpi_reduceExpectationValuesEntry(ev);
    index++;
}


void ObserverCPU::mpi_reduceExpectationValuesEntry(ExpectationValuesEntry_t &ev)
{
    unsigned numWann = sc.numWann;
    std::vector<double> sumEntries;
    std::vector<double> minEntries;
    std::vector<double> maxEntries;
    appendPackedExpectationValuesEntry(ev, sumEntries, minEntries, maxEntries);
    mpi_sum(sumEntries, sc.param.mpi);
    mpi_min(minEntries, sc.param.mpi);
    mpi_max(maxEntries, sc.param.mpi);
    unsigned sumI = 0, minI = 0, maxI = 0;
    unpackExpectationValuesEntry(ev, sumEntries, sumI, minEntries, minI, maxEntries, maxI);
}

void ObserverCPU::mpi_reduceExpectationValues()
{
    unsigned numWann = sc.numWann;
    std::vector<double> sumEntries;
    std::vector<double> minEntries;
    std::vector<double> maxEntries;
    for(unsigned evi=0; evi<expValues.size(); evi++)
        appendPackedExpectationValuesEntry(expValues[evi], sumEntries, minEntries, maxEntries);
    mpi_sum(sumEntries, sc.param.mpi);
    mpi_min(minEntries, sc.param.mpi);
    mpi_max(maxEntries, sc.param.mpi);
    unsigned sumI = 0, minI = 0, maxI = 0;
    for(unsigned evi=0; evi<expValues.size(); evi++)
        unpackExpectationValuesEntry(expValues[evi], sumEntries, sumI, minEntries, minI, maxEntries, maxI);
}

void ObserverCPU::normalizeExpectationValues(){
    for(unsigned evi=0; evi<expValues.size(); evi++){
        ExpectationValuesEntry_t &ev = expValues[evi];
        for(unsigned rid=0; rid<region.size(); rid++){
            double rw = ev.region[rid].kPointWeightSum;
            for(unsigned i=0; i<3; i++)
                ev.region[rid].j[i] /= rw;
            for(unsigned i=0; i<ev.region[rid].occupationHam.size(); i++){
                ev.region[rid].occupationHam[i] /= rw;
                ev.region[rid].occupationWan[i] /= rw;
            }
            for(unsigned i=0; i<ev.region[rid].wan.size(); i++)
                ev.region[rid].wan[i] /= rw;
        }
    }
}
