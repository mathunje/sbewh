#ifndef SBE_WH_DIAGONALIZATOR_SINGLE_H
#define SBE_WH_DIAGONALIZATOR_SINGLE_H

#include "sbewhConfig.h"

#include<algorithm>
#include<cassert>
#include<cmath>
#include<complex>
#include<vector>


#ifdef SBE_WH_LAPACK
#include<mkl_lapacke.h>
#endif

#ifdef SBE_WH_OPENMP
#include<omp.h>
#endif

enum DiagMode { jacobi, lzheev, lzheevd };

class DiagonalizatorSingle{
    const double eps = 1e-15;
    DiagMode diagMode;
    std::vector<std::complex<double>> works;
    std::vector<double> rworks;
    std::vector<int> iworks;
    std::vector<double> eigVals;
    unsigned matSize;

    unsigned sweepCount = 10;
    inline std::pair<double, double> csHalf(double y, double x);
    inline void jacobiRotate(unsigned p, unsigned q, std::complex<double> *A, std::complex<double> *U);
    inline void diagonalizeMatrix(int matId);

    void initJacobi();
    int diagJacobi(std::complex<double> * A, double * eigVals);
#ifdef SBE_WH_LAPACK
    void initZheev();
    int diagZheev(std::complex<double> * A, double * eigVals);
    void initZheevd();
    int diagZheevd(std::complex<double> * A, double * eiVals);
#endif
public:
    DiagonalizatorSingle(DiagMode dm, unsigned matSize);
    void setSweepCount(int sc) { sweepCount = sc; };
    int diag(std::complex<double> * A, double * eigenValues);
    double getError(const std::complex<double> *U, const std::complex<double> * origA) const; // (sum offdiag / sum diag) of uAu^\dagger
};

#endif
