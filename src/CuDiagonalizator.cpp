#include "CuDiagonalizator.h"

CuDiagonalizator::CuDiagonalizator(cudaStream_t nStream, unsigned matCount, unsigned matSize, unsigned offset)
    : stream(nStream), matCount(matCount), matSize(matSize), offset(offset)
{
    assert(stream);

    CUSOLVER_CHECK(cusolverDnCreate(&cusolverH));
    CUSOLVER_CHECK(cusolverDnSetStream(cusolverH, stream));
    CUSOLVER_CHECK(cusolverDnCreateSyevjInfo(&heevj_params));
    CUSOLVER_CHECK(cusolverDnXsyevjSetMaxSweeps(heevj_params, sweepCount));
    CUSOLVER_CHECK(cusolverDnXsyevjSetTolerance(heevj_params, tol));
    CUSOLVER_CHECK(cusolverDnXsyevjSetSortEig(heevj_params, 0)); // do not sort eigenvalues
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_eigvals), sizeof(double) * matSize * matCount));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(int) *  matCount));
    // query workspace
    CUSOLVER_CHECK(cusolverDnZheevjBatched_bufferSize(cusolverH, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                                                      matSize, nullptr, matSize, nullptr,
                                                      &cuLwork, heevj_params, matCount));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_work), sizeof(cuDoubleComplex) * cuLwork));

}

CuDiagonalizator::~CuDiagonalizator()
{
    CUDA_CHECK(cudaFree(d_eigvals));
    CUDA_CHECK(cudaFree(d_info));
    CUDA_CHECK(cudaFree(d_work));
    CUSOLVER_CHECK(cusolverDnDestroySyevjInfo(heevj_params));
    CUSOLVER_CHECK(cusolverDnDestroy(cusolverH));
}

void CuDiagonalizator::setTolerance(double newTol)
{
    tol = newTol;
    CUSOLVER_CHECK(cusolverDnXsyevjSetTolerance(heevj_params, tol));
}

void CuDiagonalizator::setSweepCount(int sc)
{
    sweepCount = sc;
    CUSOLVER_CHECK(cusolverDnXsyevjSetMaxSweeps(heevj_params, sweepCount));
}

void CuDiagonalizator::diag(cuDoubleComplex *A)
{
    CUSOLVER_CHECK(cusolverDnZheevjBatched(cusolverH, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                                           matSize, A + offset, matSize, d_eigvals,
                                           d_work, cuLwork,
                                           d_info, heevj_params, matCount));
}

void CuDiagonalizator::setSortEig(bool enable)
{
    int sortEig = enable ? 1 : 0;
    CUSOLVER_CHECK(cusolverDnXsyevjSetSortEig(heevj_params, sortEig));
}

std::vector<int> CuDiagonalizator::getDiagInfo() const
{
    std::vector<int> cuInfo(matCount, 0);
    CUDA_CHECK(cudaMemcpyAsync(cuInfo.data(), d_info, sizeof(int) * matCount,
                               cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    return cuInfo;
}

std::vector<double> CuDiagonalizator::getEigenValues(bool sorted) const
{
    std::vector<double> evs(matSize * matCount, 0);
    CUDA_CHECK(cudaMemcpyAsync(evs.data(), d_eigvals, sizeof(double) * evs.size(),
                               cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    if ( sorted )
        for(unsigned i=0; i<matCount; i++)
            std::sort(evs.data() + i * matSize, evs.data() + (i+1) * matSize);
    return evs;
}

std::pair<std::vector<double>, std::vector<double>> CuDiagonalizator::getBandLimits(std::vector<double> &evs) const
{
    std::vector<double> lower(matSize), upper(matSize);
    for(unsigned i=0; i<matSize; i++){
        lower[i] = evs[i];
        upper[i] = evs[i];
    }
    for(unsigned mc=1; mc<matCount; mc++)
        for(unsigned i=0; i<matSize; i++){
            double e = evs[mc*matSize+i];
            if ( e < lower[i] )
                lower[i] = e;
            if ( e > upper[i] )
                upper[i] = e;
        }
    return {lower, upper};
}
