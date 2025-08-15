#ifndef SBE_WH_CU_DIAGONALIZATOR_H
#define SBE_WH_CU_DIAGONALIZATOR_H

#include<cassert>

#include <cuda_runtime.h>
#include <cusolverDn.h>

#include<vector>
#include<algorithm>

#include "cuUtils.h"

class CuDiagonalizator {
    cudaStream_t stream = nullptr;
    unsigned matCount;
    unsigned matSize;
    unsigned offset;
    double tol = 1e-15;
    int sweepCount = 10;
    cusolverDnHandle_t cusolverH = nullptr;
    syevjInfo_t heevj_params = nullptr;
    double * d_eigvals = nullptr;
    int * d_info = nullptr;
    cuDoubleComplex * d_work = nullptr;
    int cuLwork;
public:
    CuDiagonalizator(cudaStream_t nStream, unsigned matCount, unsigned matSize, unsigned offset=0);
    ~CuDiagonalizator();
    void setSortEig(bool enable);
    void setSweepCount(int sc);
    inline int getSweepCount() const { return sweepCount; }
    void setTolerance(double newTol);
    inline double getTolerance() const { return tol; }
    void diag(cuDoubleComplex* A);
    std::vector<int> getDiagInfo() const;
    inline const double * getEigenValues_d() const { return d_eigvals; };
    inline const double * getEigenValues_d(unsigned matId) const { return d_eigvals + matId * matSize; }

    std::vector<double> getEigenValues(bool sorted=false) const;
    std::pair<std::vector<double>, std::vector<double>> getBandLimits(std::vector<double> &evs) const;
};


#endif
