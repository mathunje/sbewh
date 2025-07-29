#include "DiagonalizatorMultiple.h"


#ifdef SBE_WH_OPENMP
DiagonalizatorMultiple::DiagonalizatorMultiple(DiagMode dm, unsigned matCount, unsigned matSize, unsigned stride, unsigned offset) :
     matCount(matCount), matSize(matSize), stride(stride), offset(offset)
{
    for(unsigned i=0; i<omp_get_max_threads(); i++)
        diagSingle.push_back(DiagonalizatorSingle(dm, matSize));
    info.resize(omp_get_max_threads());
#else
DiagonalizatorMultiple::DiagonalizatorMultiple(DiagMode dm, unsigned matCount, unsigned matSize, unsigned stride, unsigned offset) :
    matCount(matCount), matSize(matSize), stride(stride), offset(offset), diagSingle(dm, matSize)
{
#endif
    assert ( matSize*matSize <= stride );
    As = nullptr;
    eigVals.resize(matCount * matSize);
}

int DiagonalizatorMultiple::diag(std::complex<double> *AsNew)
{
    As = AsNew;
#ifdef SBE_WH_OPENMP
    for(int i=0; i<omp_get_max_threads(); i++)
        info[i] = 0;
    #pragma omp parallel for schedule(static)
    for(unsigned i=0; i<matCount; i++){
        int tn = omp_get_thread_num();
        if ( info[tn] )
            continue;
        info[tn] = diagSingle[tn].diag(getAptr(i), eigVals.data()+i*matSize);
    }
    for(int i=0; i<omp_get_max_threads(); i++)
        if ( info[i] )
            return info[i];
#else
    for(unsigned i=0; i<matCount; i++){
        int info = diagSingle.diag(getAptr(i), eigVals.data()+i*matSize);
        if ( info )
            return info;
    }
#endif
    return 0;
}

double DiagonalizatorMultiple::getError(unsigned matId, const std::complex<double> * origA) const
{
    const std::complex<double> * U = getAptr(matId);
#ifdef SBE_WH_OPENMP
    return diagSingle[0].getError(U, origA);
#else
    return diagSingle.getError(U, origA);
#endif
}

double DiagonalizatorMultiple::getMaxError(const std::complex<double> *origA) const
{
    double maxError = 0;
    for(unsigned i=0; i<matCount; i++){
        double cError = getError(i, origA + i * stride + offset);
        if ( cError > maxError )
            maxError = cError;
    }
    return maxError;
}

std::pair<std::vector<double>, std::vector<double>> DiagonalizatorMultiple::getBandLimits() const
{
    std::vector<double> bandMin(matSize, std::numeric_limits<double>::max() ),
                        bandMax(matSize, std::numeric_limits<double>::lowest() );
    std::vector<double> temp(matSize);
    for(unsigned k=0; k<eigVals.size(); k+=matSize){
        std::copy(eigVals.begin()+k, eigVals.begin()+k+matSize, temp.begin());
        std::sort(temp.begin(), temp.end());
        for(unsigned i=0; i<matSize; i++){
            double e = temp[i];
            if ( e < bandMin[i] )
                bandMin[i] = e;
            if ( e > bandMax[i] )
                bandMax[i] = e;
        }
    }
    return {bandMin, bandMax};
}

void DiagonalizatorMultiple::setSweepCount(int sc)
{
#ifdef SBE_WH_OPENMP
    return diagSingle[0].setSweepCount(sc);
#else
    return diagSingle.setSweepCount(sc);
#endif
}
