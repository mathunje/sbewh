#include "CuObserver.h"

CuObserver::CuObserver(cudaStream_t stream, const StateContext &sc,
                           CuWhFourierTransformHkFull &wft, cuDoubleComplex * d_H, CuDiagonalizator &cud) :
    ObserverBase(sc), cud(cud), wft(wft), stream(stream), d_H(d_H), eim(sc)
{
    CUBLAS_CHECK ( cublasCreate(&cublasHandle) );
    CUBLAS_CHECK ( cublasSetStream(cublasHandle, stream) );
    unsigned numKpts = sc.param.prop.Nk[0] * sc.param.prop.Nk[1] * sc.param.prop.Nk[2];
    unsigned numWann = sc.numWann;
    unsigned opSize = numKpts * sc.numWann * sc.numWann;
    unsigned numRegions = region.size();
    CUDA_CHECK ( cudaMalloc(reinterpret_cast<void**>(&d_temp1), sizeof(cuDoubleComplex) * opSize) );
    CUDA_CHECK ( cudaMalloc(reinterpret_cast<void**>(&d_temp2), sizeof(cuDoubleComplex) * opSize) );
    numOps = eim.numOps;
    CUDA_CHECK ( cudaMalloc(reinterpret_cast<void**>(&d_weights), sizeof(double) * numKpts * numRegions ) );
    if ( sc.param.ksr.storeDensityMatrix ){
        CUDA_CHECK ( cudaMalloc(reinterpret_cast<void**>(&d_wan), sizeof(cuDoubleComplex) * numWann * numWann * numRegions ) );
        CUDA_CHECK ( cudaMalloc(reinterpret_cast<void**>(&d_complexWeights), sizeof(cuDoubleComplex) * numKpts * numRegions ) );
        auto begin = thrust::device_ptr<cuDoubleComplex>(d_complexWeights);
        auto end = thrust::device_ptr<cuDoubleComplex>(d_complexWeights + numKpts * numRegions );
        thrust::fill(begin, end, cuDoubleComplex(0, 0));
        wan.resize(numWann * numWann * numRegions);
    }
    CUDA_CHECK ( cudaMalloc(reinterpret_cast<void**>(&d_kexp), sizeof(double) * numKpts * numOps)  );
    CUDA_CHECK ( cudaMalloc(reinterpret_cast<void**>(&d_regionOps), sizeof(double) * numRegions * numOps) );
    regionOps.resize(numOps * numRegions);
    d_kBaseShifts = wft.createKshifts(sc.kGlobalShift_au);
    CUDA_CHECK ( cudaMalloc(reinterpret_cast<void**>(&d_observerRegion), sizeof(ObserverRegion_t) * numRegions) );
    CUDA_CHECK ( cudaMemcpy(d_observerRegion, region.data(), sizeof(ObserverRegion_t) * numRegions, cudaMemcpyHostToDevice) );

    specialValues_h.resize(sc.si.specialIndexCount);
}

CuObserver::~CuObserver()
{
    CUDA_CHECK ( cudaFree(d_temp1) );
    CUDA_CHECK ( cudaFree(d_temp2) );
    CUDA_CHECK ( cudaFree(d_weights) );
    if ( sc.param.ksr.storeDensityMatrix ){
        CUDA_CHECK ( cudaFree(d_wan) );
        CUDA_CHECK ( cudaFree(d_complexWeights) );
    }
    CUDA_CHECK ( cudaFree(d_kexp) );
    CUDA_CHECK ( cudaFree(d_regionOps) );
    CUDA_CHECK ( cudaFree(d_kBaseShifts) );
    CUDA_CHECK ( cudaFree(d_observerRegion) );
    CUBLAS_CHECK ( cublasDestroy(cublasHandle) );
}


__global__ void fillOnes(double *dst, unsigned size)
{
    unsigned id = threadIdx.x + blockDim.x * blockIdx.x;
    if ( id >= size)
        return;
    dst[id] = 1;
}

__global__ void calcWeights(double *weights, ObserverRegion_t * region, double * kBaseShifts, unsigned kPaddedDim, unsigned numKpts,
                            double A0, double A1, double A2)
{
    double A[3] = {A0, A1, A2};
    unsigned id = threadIdx.x + blockDim.x * blockIdx.x;
    if ( id >= numKpts)
        return;
    unsigned rid = blockIdx.y;
    double * k = kBaseShifts + kPaddedDim * id;
    double w = 1;
    for(unsigned dir=0; dir<3; dir++){
        double kRef = (region[rid].movingFrame ? k[dir] : k[dir] - A[dir] );
        kRef -= region[rid].k[dir];
        double cPos = kRef - std::round(kRef);
        w *= std::exp( -cPos * cPos * region[rid].invSigma2);
    }
    weights[id * gridDim.y + rid] = w;
}

__global__ void calcJk(double *dst, cuDoubleComplex *T, const cuDoubleComplex *rho, const cuDoubleComplex *D, const cuDoubleComplex * dHdK,
                       unsigned numWann, unsigned maxDim, unsigned numKpts)
{
    unsigned id = threadIdx.x + blockDim.x * blockIdx.x;
    if ( id >= numKpts)
        return;
    unsigned baseSize = numKpts * numWann * numWann;
    for(unsigned dir=0; dir<maxDim; dir++){
        const cuDoubleComplex * Dk = D + dir * baseSize + id * numWann * numWann;
        const cuDoubleComplex * dHk = dHdK + dir * baseSize + id * numWann * numWann;
        const cuDoubleComplex * rhok = rho + id * numWann * numWann;
        const cuDoubleComplex * Tk = T + id * numWann * numWann;
        double sum = 0.0;
        for(unsigned a=0; a<numWann; a++)
            for(unsigned b=0; b<numWann; b++){
                sum += Dk[a*numWann+b].x * Tk[a*numWann+b].x - Dk[a*numWann+b].y * Tk[a*numWann+b].y; // real part of trace(DT^T)
                sum -= rhok[a*numWann+b].x * dHk[b*numWann+a].x - rhok[a*numWann+b].y * dHk[b*numWann+a].y; // real part of trace( dHdK rho)
            }
        dst[id+numKpts*dir] = sum;
    }
}

__global__ void extractOccWan(double *dst, const cuDoubleComplex * rhoW, unsigned numWann, unsigned numKpts)
{
    unsigned id = threadIdx.x + blockDim.x * blockIdx.x;
    if ( id >= numKpts)
        return;
    for(unsigned i=0; i<numWann; i++)
        dst[i * numKpts + id] = rhoW[id * numWann * numWann + i * (numWann+1)].x;
}

__global__ void extractOccHam(double *dst, cuDoubleComplex * rhoH, unsigned numWann, unsigned numKpts)
{
    unsigned id = threadIdx.x + blockDim.x * blockIdx.x;
    if ( id >= numKpts)
        return;
    for(unsigned i=0; i<numWann; i++)
        dst[i + id * numWann] = rhoH[id * numWann * numWann + i * (numWann+1)].x;
}

__global__ void smearOcc(double *smearedOcc, const double *occ, const double * energy, double occupationSmearingWidth, unsigned numWann, unsigned numKpts)
{
    unsigned id = threadIdx.x + blockDim.x * blockIdx.x;
    if ( id >= numKpts * numWann)
        return;
    unsigned a = id % numWann;
    unsigned kpt = id / numWann;
    if ( occupationSmearingWidth < 0 ){
        smearedOcc[kpt + a*numKpts] = occ[a + kpt*numWann];
    } else {
        double sumW = 0;
        double res = 0;
        double selfE = energy[a + kpt * numWann];
        for(unsigned b=0; b<numWann; b++){
            double dE = (selfE - energy[b + kpt*numWann]) / occupationSmearingWidth;
            double w = exp(-dE*dE);
            res += w * occ[b + kpt*numWann];
            sumW +=w;
        }
        smearedOcc[kpt + a*numKpts] = res / sumW;
    }
}


void CuObserver::operator()(const GPUstate &state, double t)
{
    /*
     * NOTE: As the cublas functions assume column-major ordering, contary to CPU implementation,
     *       the matrix operations are transposed and performed in reversed order compared to CPU code.
     */
    basicLog(t);
    unsigned numWann = sc.numWann;
    unsigned maxDim = sc.dim;
    unsigned numKpts = sc.param.prop.Nk[0] * sc.param.prop.Nk[1] * sc.param.prop.Nk[2];
    unsigned numRegions = region.size();
    unsigned baseCount = numWann * numWann * numKpts;
    bool storeRegionDensity = getStoreRegionDensity( expValues.size() );
    ExpectationValuesEntry_t ev = prepareExpValue(t, storeRegionDensity);

    const cuDoubleComplex * statePtr = reinterpret_cast<const cuDoubleComplex*>(thrust::raw_pointer_cast(state.data()));
    const cuDoubleComplex * rho = statePtr + sc.si.specialIndexCount;

    CUDA_CHECK(cudaMemcpyAsync(specialValues_h.data(), statePtr, sizeof(std::complex<double>) * sc.si.specialIndexCount,
                                 cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));

    GeomVector3d A;
    for(unsigned i=0; i<3; i++){
        A[i] = std::real(specialValues_h[sc.si.A[i]]);
        ev.A[i] = std::real(A[i]);
    }

    GeomVector3d kOff = sc.kGlobalShift_au + A;
    wft.transform(d_H, kOff);
    CUDA_CHECK(cudaStreamSynchronize(stream));

    cuDoubleComplex * d_H0 = d_H + TB_FULL_OP_POS(H0) * baseCount;
    // temp1 = -i * [H_0, rho]^T = -i * [rho^T, H_0^T]
    cuDoubleComplex alpha = {.x=0, .y=-1}, beta = {.x=0, .y=0};
    CUBLAS_CHECK (
            cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_T, CUBLAS_OP_T,
                  numWann, numWann, numWann,
                  &alpha,
                  rho, numWann, numWann*numWann,
                  d_H0, numWann, numWann*numWann,
                  &beta,
                  d_temp1, numWann, numWann*numWann,
                  numKpts)
                 );
    cuDoubleComplex alpha2 = {.x=0, .y=1};
    cuDoubleComplex beta2 = {.x=1, .y=0};
    CUBLAS_CHECK ( cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_T, CUBLAS_OP_T,
                  numWann, numWann, numWann,
                  &alpha2,
                  d_H0, numWann, numWann*numWann,
                  rho, numWann, numWann*numWann,
                  &beta2,
                  d_temp1, numWann, numWann*numWann,
                  numKpts)
                 );
    unsigned nThreadsPerBlock = 256;
    dim3 nBlocks( (numKpts+nThreadsPerBlock-1) / nThreadsPerBlock,
                  numRegions,
                  1
                 );
    GeomVector3d Afrac = wft.toFracK(A);
    calcWeights<<< nBlocks, nThreadsPerBlock, 0, stream >>>
            (d_weights, d_observerRegion, d_kBaseShifts, wft.getKpaddedDim(), numKpts,
             Afrac[0], Afrac[1], Afrac[2]);

    fillOnes<<<(numKpts + nThreadsPerBlock-1)/nThreadsPerBlock, nThreadsPerBlock, 0, stream >>>
            (d_kexp + numKpts * eim.weights, numKpts);
    cuDoubleComplex * d_D = d_H + TB_FULL_OP_POS(D) * baseCount;
    cuDoubleComplex * d_dHdK = d_H + TB_FULL_OP_POS(dH0dk) * baseCount;
    calcJk<<<(numKpts + nThreadsPerBlock-1)/nThreadsPerBlock, nThreadsPerBlock, 0, stream >>>
            (d_kexp + numKpts * eim.j, d_temp1, rho, d_D, d_dHdK, numWann, maxDim, numKpts);

    cud.diag(d_H0);
    cuDoubleComplex * U = d_H0; // overriden by cud.diag
    CUDA_CHECK(cudaStreamSynchronize(stream));
    extractOccWan<<<(numKpts + nThreadsPerBlock-1)/nThreadsPerBlock, nThreadsPerBlock, 0, stream >>>
            (d_kexp + numKpts * eim.occWan, rho, numWann, numKpts);
    // d_temp2  = \rho^H  = U^\dagger \rho U
    alpha = { .x = 1.0, .y = 0.0};
    beta =  { .x = 0, .y = 0 };
    CUBLAS_CHECK (
            cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_C, CUBLAS_OP_N,
                  numWann, numWann, numWann,
                  &alpha,
                  U, numWann, numWann * numWann,
                  rho, numWann, numWann*numWann,
                  &beta,
                  d_temp1, numWann, numWann*numWann,
                  numKpts)
                 );
    CUBLAS_CHECK ( cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N,
                  numWann, numWann, numWann,
                  &alpha,
                  d_temp1, numWann, numWann * numWann,
                  U, numWann, numWann*numWann,
                  &beta,
                  d_temp2, numWann, numWann*numWann,
                  numKpts)
                 );

    double * d_temp1Occs = reinterpret_cast<double*>(d_temp1);

    extractOccHam<<<(numKpts + nThreadsPerBlock-1)/nThreadsPerBlock, nThreadsPerBlock, 0, stream >>>
            (d_temp1Occs, d_temp2, numWann, numKpts);
    smearOcc<<<(numKpts * numWann + nThreadsPerBlock-1)/nThreadsPerBlock, nThreadsPerBlock, 0, stream >>>
            (d_kexp + numKpts * eim.occHam, d_temp1Occs, cud.getEigenValues_d(), sc.param.prop.occupationSmearingWidth, numWann, numKpts);

    CUDA_CHECK(cudaStreamSynchronize(stream));
    for(unsigned i=0; i<numWann; i++){
        auto begin = thrust::device_ptr<double>(d_kexp + numKpts * (eim.occHam+i) );
        auto end = thrust::device_ptr<double>(d_kexp + numKpts * (eim.occHam+i+1) );
        ev.occupationHamMin[i] = thrust::reduce(begin, end,
                                                std::numeric_limits<double>::max(), thrust::minimum<double>());
        ev.occupationHamMax[i] = thrust::reduce(begin, end,
                                                std::numeric_limits<double>::lowest(), thrust::maximum<double>());
    }

    double alphaD = 1, betaD = 0;
    CUBLAS_CHECK (
            cublasDgemm(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N,
                  numRegions, numOps, numKpts,
                  &alphaD,
                  d_weights, numRegions,
                  d_kexp, numKpts,
                  &betaD,
                  d_regionOps, numRegions)
                 );
    CUDA_CHECK( cudaMemcpyAsync(regionOps.data(), d_regionOps, sizeof(double) * numOps * numRegions,
                                 cudaMemcpyDeviceToHost, stream));

    if ( storeRegionDensity ){
        cuDoubleComplex alphaZ {1.0, 0.0}, betaZ {0.0, 0.0};
        CUBLAS_CHECK ( cublasDcopy( cublasHandle, numKpts * numRegions, d_weights, 1, (double*)d_complexWeights, 2) );
        CUBLAS_CHECK (
                cublasZgemm(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_T,
                            numRegions, numWann*numWann, numKpts,
                            &alphaZ,
                            d_complexWeights, numRegions,
                            rho, numWann*numWann,
                            &betaZ,
                            d_wan, numRegions)
                );
        CUDA_CHECK( cudaMemcpyAsync(wan.data(), d_wan, sizeof(cuDoubleComplex) * numWann * numWann * numRegions,
                                 cudaMemcpyDeviceToHost, stream));
    }

    CUDA_CHECK(cudaStreamSynchronize(stream));

    for(unsigned rid=0; rid<numRegions; rid++){
        double w = regionOps[eim.weights * numRegions + rid];
        ev.region[rid].kPointWeightSum = w;
        for(unsigned dir=0; dir<maxDim; dir++){
            ev.region[rid].j[dir] = regionOps[(eim.j + dir) * numRegions + rid] / w;
        }
        for(unsigned bi=0; bi<numWann; bi++){
            ev.region[rid].occupationWan[bi] = regionOps[(eim.occWan + bi) * numRegions + rid] / w;
            ev.region[rid].occupationHam[bi] = regionOps[(eim.occHam + bi) * numRegions + rid] / w;
        }
        if ( storeRegionDensity ){
            for(unsigned i=0; i<numWann*numWann; i++)
                ev.region[rid].wan[i] = wan[i*numRegions + rid] / w;
        }
    }
    expValues.push_back(ev);
}
