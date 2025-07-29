#include "CuPropagator.h"

CuPropagator::CuPropagator(cudaStream_t stream, const StateContext &sc,
                           CuWhFourierTransformHk &wft, cuDoubleComplex * d_H, CuDiagonalizator &cud) :
    PropagatorBase(sc), stream(stream), wft(wft), d_H(d_H), cud(cud)
{
    CUBLAS_CHECK ( cublasCreate(&cublasHandle) );
    CUBLAS_CHECK ( cublasSetStream(cublasHandle, stream) );
    if ( sc.param.prop.relaxationF2 > 0){
        unsigned size = sc.param.prop.Nk[0] * sc.param.prop.Nk[1] * sc.param.prop.Nk[2] * sc.numWann * sc.numWann;
        CUDA_CHECK ( cudaMalloc(reinterpret_cast<void**>(&d_intermediate), sizeof(cuDoubleComplex) * size) );
        CUDA_CHECK ( cudaMalloc(reinterpret_cast<void**>(&d_rhoH), sizeof(cuDoubleComplex) * size) );
    }
    specialValues_h.resize(sc.si.specialIndexCount);
    specialValues_dsdt_h.resize(sc.si.specialIndexCount);
}


CuPropagator::~CuPropagator()
{
    if ( sc.param.prop.relaxationF2 > 0){
        CUDA_CHECK ( cudaFree(d_intermediate) );
        CUDA_CHECK ( cudaFree(d_rhoH) );
    }
    CUBLAS_CHECK ( cublasDestroy(cublasHandle) );
}

__global__ void fillOccupancies(cuDoubleComplex *state, const double * eigVals, double fermiLevel,
                                double gsTemp, unsigned numWann, unsigned size)
{
    unsigned id = threadIdx.x + blockDim.x * blockIdx.x;
    if ( id >= size)
        return;
    state[(id/numWann) *numWann*numWann + (id%numWann)*(numWann+1)] =
        cuDoubleComplex(1 / ( 1 + exp( (eigVals[id]-fermiLevel)/gsTemp ) ), 0);
}


GPUstate CuPropagator::calcInitialState(double fermiLevel)
{
    unsigned nKpts = sc.param.prop.Nk[0] * sc.param.prop.Nk[1] * sc.param.prop.Nk[2];
    unsigned numWann = sc.numWann;
    wft.transform(d_H, sc.kGlobalShift_au, {0, 0, 0});
    unsigned baseSize = nKpts * numWann*numWann;
    cuDoubleComplex * d_H0 = d_H + TB_OP_POS(H0) * baseSize;
    cuDoubleComplex * d_temp = d_H + TB_OP_POS(HE) * baseSize;
    cud.diag(d_H0);
    cuDoubleComplex * U = d_H0; // overriden by cud.diag
    // 2 doubles form one complex number
    GPUstate state( 2 * (sc.si.specialIndexCount + baseSize) , 0); // sets A(t_start) = 0, too.
    cuDoubleComplex * statePtr = reinterpret_cast<cuDoubleComplex*>(thrust::raw_pointer_cast(state.data()))
                                 + sc.si.specialIndexCount;
    unsigned nThreadsPerBlock = 256;
    fillOccupancies<<< (baseSize + nThreadsPerBlock-1)/nThreadsPerBlock, nThreadsPerBlock, 0, stream >>>
        (statePtr, cud.getEigenValues_d(), fermiLevel, sc.param.tb.gsTemp, numWann, numWann * nKpts);
    cuDoubleComplex alpha = { .x = 1.0, .y = 0.0}, beta =  { .x = 0, .y = 0 };
    CUBLAS_CHECK (
            cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N,
                  numWann, numWann, numWann,
                  &alpha,
                  U, numWann, numWann * numWann,
                  statePtr, numWann, numWann*numWann,
                  &beta,
                  d_temp, numWann, numWann*numWann,
                  nKpts)
                 );
    CUBLAS_CHECK ( cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_C,
                  numWann, numWann, numWann,
                  &alpha,
                  d_temp, numWann, numWann * numWann,
                  U, numWann, numWann*numWann,
                  &beta,
                  statePtr, numWann, numWann*numWann,
                  nKpts)
                 );
    CUDA_CHECK(cudaStreamSynchronize(stream));
    return state;
}

__global__ void toDephasing(cuDoubleComplex *state, const double * energy, const double f2,
                            const double soothingWidth, unsigned numWann, unsigned matCount)
{
    unsigned id = threadIdx.x + blockDim.x * blockIdx.x;
    unsigned matId = id / numWann / numWann;
    if ( matId >= matCount )
        return;
    unsigned a = (id / numWann) % numWann;
    unsigned b = id % numWann;
    if ( soothingWidth < 0 ){
        double preFact = - (a != b) * f2;
        state[id] = cuDoubleComplex(preFact * state[id].x, preFact * state[id].y);
    } else {
        double dE = (energy[matId*numWann+a] - energy[matId*numWann+b]) / soothingWidth;
        double preFact = - f2* (1 - exp(-dE*dE));
        state[id] = cuDoubleComplex(preFact * state[id].x, preFact * state[id].y);
    }
}


void CuPropagator::operator() (const GPUstate &state, GPUstate &dsdt, const double t)

{
    /*
     * NOTE: As the cublas functions assume column-major ordering, contary to CPU implementation,
     *       the matrix operations are transposed and performed in reversed order compared to CPU code.
     */
    auto startTime = std::chrono::steady_clock::now();
    unsigned numWann = sc.numWann;
    unsigned nKpts = sc.param.prop.Nk[0] * sc.param.prop.Nk[1] * sc.param.prop.Nk[2];
    unsigned baseSize = nKpts * numWann * numWann;
    const cuDoubleComplex * statePtr = reinterpret_cast<const cuDoubleComplex*>(thrust::raw_pointer_cast(state.data()));
    const cuDoubleComplex * rho = statePtr + sc.si.specialIndexCount;
    cuDoubleComplex * dsdtPtr = reinterpret_cast<cuDoubleComplex*>(thrust::raw_pointer_cast(dsdt.data()));
    cuDoubleComplex * dRho = dsdtPtr + sc.si.specialIndexCount;
    // thrust::fill(dsdt.begin(), dsdt.end(), 0);
    GeomVector3d E = sc.getE(t);
    CUDA_CHECK(cudaMemcpyAsync(specialValues_h.data(), statePtr, sizeof(std::complex<double>) * sc.si.specialIndexCount,
                                 cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));

    GeomVector3d A({ 0, 0, 0 });
    for(unsigned i=0; i<sc.dim; i++){
        A[i] = std::real(specialValues_h[sc.si.A[i]]);
        specialValues_dsdt_h[sc.si.A[i]] = -E[i];
    }
    CUDA_CHECK(cudaMemcpyAsync(dsdtPtr, specialValues_dsdt_h.data(), sizeof(std::complex<double>) * sc.si.specialIndexCount,
                                 cudaMemcpyHostToDevice, stream));

    GeomVector3d kOff = sc.kGlobalShift_au + A;
    wft.transform(d_H, kOff, E);
    CUDA_CHECK(cudaStreamSynchronize(stream));


    cuDoubleComplex * d_HE = d_H + TB_OP_POS(HE) * baseSize;
    // coherent propagation, i.e. i[H, rho]
    cuDoubleComplex alpha = {.x=0, .y=1}, beta = {.x=0, .y=0};
    CUBLAS_CHECK (
            cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N,
                  numWann, numWann, numWann,
                  &alpha,
                  d_HE, numWann, numWann*numWann,
                  rho, numWann, numWann*numWann,
                  &beta,
                  dRho, numWann, numWann*numWann,
                  nKpts)
                 );
    cuDoubleComplex alpha2 = {.x=0, .y=-1};
    cuDoubleComplex beta2 = {.x=1, .y=0};
    CUBLAS_CHECK ( cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N,
                  numWann, numWann, numWann,
                  &alpha2,
                  rho, numWann, numWann*numWann,
                  d_HE, numWann, numWann*numWann,
                  &beta2,
                  dRho, numWann, numWann*numWann,
                  nKpts)
                 );

    // dephasing in Hamiltonian gauge
    if ( sc.param.prop.relaxationF2 > 0){
        cuDoubleComplex * d_H0 = d_H + TB_OP_POS(H0) * baseSize;
        cud.diag(d_H0);
        CUDA_CHECK(cudaStreamSynchronize(stream));
        cuDoubleComplex * U = d_H0; // overriden by cud.diag
        alpha = { .x = 1.0, .y = 0.0};
        beta = { .x = 0, .y = 0 };
        CUBLAS_CHECK (
                cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_C, CUBLAS_OP_N,
                      numWann, numWann, numWann,
                      &alpha,
                      U, numWann, numWann * numWann,
                      rho, numWann, numWann*numWann,
                      &beta,
                      d_intermediate, numWann, numWann*numWann,
                      nKpts)
                     );
        CUBLAS_CHECK ( cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N,
                      numWann, numWann, numWann,
                      &alpha,
                      d_intermediate, numWann, numWann * numWann,
                      U, numWann, numWann*numWann,
                      &beta,
                      d_rhoH, numWann, numWann*numWann,
                      nKpts)
                     );
        // d_rhoH -> dephasing in hamiltonian gauge
        unsigned nThreadsPerBlock = 256;
        toDephasing<<< (baseSize + nThreadsPerBlock-1)/nThreadsPerBlock, nThreadsPerBlock, 0, stream >>>
            (d_rhoH, cud.getEigenValues_d(), sc.param.prop.relaxationF2, sc.param.prop.soothingWidth, numWann, nKpts);
        CUDA_CHECK(cudaStreamSynchronize(stream));
        // transform to Wannier gauge and add to dRho.
        CUBLAS_CHECK (
                cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N,
                      numWann, numWann, numWann,
                      &alpha,
                      U, numWann, numWann * numWann,
                      d_rhoH, numWann, numWann*numWann,
                      &beta,
                      d_intermediate, numWann, numWann*numWann,
                      nKpts)
                     );
        cuDoubleComplex betaAdd = {.x = 1, .y = 0};
        CUBLAS_CHECK ( cublasZgemmStridedBatched(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_C,
                      numWann, numWann, numWann,
                      &alpha,
                      d_intermediate, numWann, numWann * numWann,
                      U, numWann, numWann*numWann,
                      &betaAdd,
                      dRho, numWann, numWann*numWann,
                      nKpts)
                     );
    }
    CUDA_CHECK(cudaStreamSynchronize(stream));
    if ( sc.param.general.timeIt )
        timeAndPrint(startTime);
}
