#ifndef SBE_WH_DIAGONALIZATOR_MULTIPLE_H
#define SBE_WH_DIAGONALIZATOR_MULTIPLE_H

#include "sbewhConfig.h"

#include<algorithm>
#include<cassert>
#include<cmath>
#include<complex>
#include<vector>


#include "DiagonalizatorSingle.h"

#ifdef SBE_WH_OPENMP
#include<omp.h>
#endif

/*
 * Diagonalizes multiple matrices at once.
 * - assumes COL_MAJOR ordering.
 */

class DiagonalizatorMultiple {
    std::complex<double> * As;
    unsigned matCount;
    unsigned matSize;
    unsigned stride;
    unsigned offset;
#ifdef SBE_WH_OPENMP
    std::vector<int> info;
    std::vector<DiagonalizatorSingle> diagSingle;
#else
    DiagonalizatorSingle diagSingle;
#endif
    std::vector<double> eigVals;
    std::complex<double> * getAptr(unsigned matId) { return As + stride * matId + offset; }
    const std::complex<double> * getAptr(unsigned matId) const { return As + stride * matId + offset; }

public:
    DiagonalizatorMultiple(DiagMode dm, unsigned matCount, unsigned matSize, unsigned stride, unsigned offset=0);
    void setSweepCount(int sc);
    int diag(std::complex<double> * AsNew);
    // refers to last diagonalized As
    double getError(unsigned matId, const std::complex<double> * origA) const; // (sum offdiag / sum diag) of uAu^\dagger
    double getMaxError(const std::complex<double> *origA) const;
    inline const std::vector<double> & getEigenValues() const { return eigVals; };
    inline const double * getEigenValues(unsigned matId) const { return eigVals.data() + matId * matSize; }
    std::pair<std::vector<double>, std::vector<double>> getBandLimits() const;
};

#endif
