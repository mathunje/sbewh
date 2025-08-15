#include "CuWhFourierTransform.h"

CuWhFourierTransformBase::CuWhFourierTransformBase(cudaStream_t stream,
                    const TightBindingParameter_t &tb, const FourierTransformParameter_t &ftp,
                    const std::array<unsigned, 3> &Nr, unsigned numOps) :
    stream(stream), tb(tb), ftp(ftp), Nf(ftp.Nf), Nr(Nr), numOps(numOps), minCellIndices(tb.minCellIndices),
    numWann(tb.numWann)
{
    const auto &Rl = tb.latticeVectors;
    CUDA_CHECK( cudaMalloc(reinterpret_cast<void**>(&d_meshShifts), sizeof(double) * kPaddedDim * Nr[0]*Nr[1]*Nr[2]) );
    CUFFT_CHECK( cufftCreate(&plan) );
    int batch_size = Nr[0]*Nr[1]*Nr[2] * numWann*numWann * numOps;
    std::array<int, 3> Nfi({(int)Nf[0], (int)Nf[1], (int)Nf[2]});
    CUFFT_CHECK( cufftPlanMany(&plan, Nfi.size(), Nfi.data(),
                              nullptr, 1, 0, // *inembed, istride, idist
                              nullptr, 1, 0, // *onembed, ostride, odist
                              CUFFT_Z2Z, batch_size) );
    CUDA_CHECK( cudaMalloc(reinterpret_cast<void**>(&d_fftData), sizeof(cuDoubleComplex)* batch_size * Nf[0]*Nf[1]*Nf[2] ) );
    CUFFT_CHECK( cufftSetStream(plan, stream) );

    unsigned fftSize = Nf[0] * Nf[1] * Nf[2];
    std::vector<std::complex<double>> rOps( fftSize * numWann * numWann * TB_FULL_OP_CNT);
    for(const auto [ci, Hmat] : tb.realSpaceOps){
        CellIndex gridCi = ci - tb.minCellIndices;
        unsigned fftIndex = gridCi[2] + Nf[2] * ( gridCi[1] + Nf[1] * gridCi[0]);
        for(unsigned m=0; m<numWann; m++)
            for(unsigned n=0; n<numWann; n++)
                rOps[(( TB_FULL_OP_POS(H0) * numWann + m ) * numWann + n) * fftSize + fftIndex] = Hmat.H(m, n);
        GeomVector3d pos = (double)ci.at(0) * Rl[0] + (double)ci.at(1) * Rl[1] + (double)ci.at(2) * Rl[2];
        for(unsigned d=0; d<3; d++)
            for(unsigned m=0; m<tb.numWann; m++)
                for(unsigned n=0; n<tb.numWann; n++){
                    rOps[(( (TB_FULL_OP_POS(D)+d)  * numWann + m ) * numWann + n) * fftSize + fftIndex] = Hmat.D[d](m, n);
                    rOps[(( (TB_FULL_OP_POS(dH0dk)+d)  * numWann + m ) * numWann + n) * fftSize + fftIndex] = std::complex<double>{0, 1} * pos[d]  * Hmat.H(m, n);
                }
    }
    CUDA_CHECK( cudaMalloc(reinterpret_cast<void**>(&d_realSpaceGridOps), sizeof(std::complex<double>)*rOps.size()) );
    CUDA_CHECK( cudaMemcpyAsync(d_realSpaceGridOps, rOps.data(), sizeof(std::complex<double>) * rOps.size(),
                                 cudaMemcpyHostToDevice, stream) );
}

std::vector<GeomVector3d> CuWhFourierTransformBase::getKmesh(const GeomVector3d &kOffset_au) const
{
    const auto &Rl = tb.latticeVectors;
    double preFact = 1  / ( Rl[0].dot(Rl[1].cross(Rl[2])));
    GeomVector3d kBase[3] = { preFact * Rl[1].cross(Rl[2]),
                              preFact * Rl[2].cross(Rl[0]),
                              preFact * Rl[0].cross(Rl[1])
                            };
    GeomVector3d kFracOffset = toFracK(kOffset_au);
    std::vector<GeomVector3d> kMesh;
    std::array<unsigned, 3> Nk({Nr[0]*Nf[0], Nr[1]*Nf[1], Nr[2]*Nf[2]});
    kMesh.reserve(Nk[0] * Nk[1] * Nk[2]);
    for(unsigned a=0; a<Nk[0]; a++)
        for(unsigned b=0; b<Nk[1]; b++)
            for(unsigned c=0; c<Nk[2]; c++)
                kMesh.push_back( a /(double)Nk[0] * kBase[0] + b/(double)Nk[1] * kBase[1] + c/(double)Nk[2] * kBase[2]
                                 + kFracOffset );
    return kMesh;
}

CuWhFourierTransformBase::~CuWhFourierTransformBase()
{
    CUDA_CHECK( cudaFree(d_meshShifts) );
    CUDA_CHECK( cudaFree(d_realSpaceGridOps) );
    CUDA_CHECK( cudaFree(d_fftData) );
    CUFFT_CHECK ( cufftDestroy(plan) );
}


__global__ void cuUpdateMeshShifts(double * mShifts, int Nr0, int Nr1, int Nr2,
                                int Nk0, int Nk1, int Nk2, double kOff0, double kOff1, double kOff2,
                                int kPaddedDim)
{
    int kIndex = blockIdx.x * blockDim.x + threadIdx.x;
    int NrSize = Nr0 * Nr1 * Nr2;
    if ( kIndex >= NrSize )
        return;
    int a = kIndex / (Nr1 * Nr2);
    int b = (kIndex / Nr2) % Nr1;
    int c = kIndex % Nr2;
    mShifts[kPaddedDim * kIndex + 0] = a / double(Nk0) + kOff0;
    mShifts[kPaddedDim * kIndex + 1] = b / double(Nk1) + kOff1;
    mShifts[kPaddedDim * kIndex + 2] = c / double(Nk2) + kOff2;
}

void CuWhFourierTransformBase::updateMeshShifts(const GeomVector3d &kFracOffset)
{
    int nThreadsPerBlock = 256;
    cuUpdateMeshShifts<<<(Nr[0]*Nr[1]*Nr[2] + nThreadsPerBlock-1)/nThreadsPerBlock, nThreadsPerBlock, 0, stream >>>
            (d_meshShifts, Nr[0], Nr[1], Nr[2], Nr[0]*Nf[0], Nr[1]*Nf[1], Nr[2]*Nf[2],
             kFracOffset.at(0), kFracOffset.at(1), kFracOffset.at(2),
             kPaddedDim);
}

__global__ void expandWithKphase(const cuDoubleComplex * rOps, cuDoubleComplex * expandedR,
                                 const double * kShifts,
                                 int Nf0, int Nf1, int Nf2,
                                 int kPaddedDim,
                                 int kOffSize)
{
    int kOffIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if ( kOffIndex >= kOffSize )
        return;
    int Rindex = blockIdx.y;
    int outerIndex = blockIdx.z;
    const double * kOff = kShifts + kPaddedDim * kOffIndex;
    int a = Rindex / (Nf1 * Nf2);
    int b = (Rindex / Nf2) % Nf1;
    int c = Rindex % Nf2;
    double phase = 2 * std::numbers::pi * ( a * kOff[0] + b * kOff[1] + c * kOff[2] );
    cuDoubleComplex eikr(cos(phase), sin(phase));
    expandedR[ Rindex + gridDim.y * ( kOffIndex + kOffSize * outerIndex ) ] =
         cuCmul(eikr, rOps[ Rindex + outerIndex * gridDim.y ]);
}

__global__ void reorderWithPhase(const cuDoubleComplex * fftData,
                                 const double * kShifts,
                                 cuDoubleComplex * out,
                                 int mci0, int mci1, int mci2,
                                 int Nf0, int Nf1, int Nf2,
                                 int Nr0, int Nr1, int Nr2,
                                 int kPaddedDim)
{
    int NfSize = Nf0 * Nf1 * Nf2;
    int NrSize = Nr0 * Nr1 * Nr2;
    int kIndex = blockIdx.x * blockDim.x + threadIdx.x;
    int kShiftIndex = kIndex / NfSize;
    if ( kShiftIndex >= NrSize )
        return;
    int kfIndex = kIndex - kShiftIndex * NfSize;
    int mnIndex = blockIdx.y;
    int opIndex = blockIdx.z;
    int a = kfIndex / (Nf1 * Nf2);
    int b = (kfIndex / Nf2) % Nf1;
    int c = kfIndex % Nf2;
    const double * kOff = kShifts + kPaddedDim * kShiftIndex;
    double phase = 2 * std::numbers::pi * (
                      mci0 * kOff[0] + mci1 * kOff[1] + mci2 * kOff[2] +
                      (a * mci0) / (double)(Nf0) +
                      (b * mci1) / (double)(Nf1) +
                      (c * mci2) / (double)(Nf2) );
    cuDoubleComplex eikr(cos(phase), sin(phase));
    out[mnIndex + gridDim.y * ( kIndex + NfSize*NrSize * opIndex)] =
            cuCmul(eikr, fftData[ kIndex + (NfSize*NrSize) * ( mnIndex + gridDim.y * opIndex)]);
}

void CuWhFourierTransformBase::applyTransform(const cuDoubleComplex * rOps, cuDoubleComplex * kOps)
{
    int Nrep = Nr[0] * Nr[1] * Nr[2];
    int Nfft = Nf[0] * Nf[1] * Nf[2];
    int nThreadsPerBlock = 256;
    dim3 nBlocks( (Nrep+nThreadsPerBlock-1) / nThreadsPerBlock,
                  Nfft,
                  numOps * numWann*numWann
                 );
    expandWithKphase<<<nBlocks, nThreadsPerBlock, 0, stream>>>
        (rOps, d_fftData, d_meshShifts,
         Nf[0], Nf[1], Nf[2],
         kPaddedDim, Nrep);
    CUFFT_CHECK(cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(d_fftData),
                                   reinterpret_cast<cufftDoubleComplex*>(d_fftData),
                                   CUFFT_INVERSE));
    nBlocks = dim3( (Nfft*Nrep + nThreadsPerBlock-1) / nThreadsPerBlock,
                    numWann * numWann,
                    numOps);
    reorderWithPhase<<<nBlocks, nThreadsPerBlock, 0, stream>>>
        (d_fftData, d_meshShifts, kOps,
         tb.minCellIndices.at(0), tb.minCellIndices.at(1), tb.minCellIndices.at(2),
         Nf[0], Nf[1], Nf[2],
         Nr[0], Nr[1], Nr[2],
         kPaddedDim);
}

cuDoubleComplex * CuWhFourierTransformBase::createDataVec() const
{
    cuDoubleComplex * d_data = nullptr;
    unsigned size = numOps * Nr[0]*Nr[1]*Nr[2] * Nf[0]*Nf[1]*Nf[2] * numWann * numWann;
    CUDA_CHECK ( cudaMalloc(reinterpret_cast<void **>(&d_data), sizeof(cuDoubleComplex) * size) );
    return d_data;
}

__global__ void populateKshifts(double * mShifts,
                                int Nr0, int Nr1, int Nr2,
                                int Nf0, int Nf1, int Nf2,
                                int kPaddedDim,
                                double kgs0, double kgs1, double kgs2)
{
    int mIndex = blockIdx.x * blockDim.x + threadIdx.x;
    int NrSize = Nr0 * Nr1 * Nr2;
    if ( mIndex >= NrSize )
        return;
    int ma = mIndex / (Nr1 * Nr2);
    int mb = (mIndex / Nr2) % Nr1;
    int mc = mIndex % Nr2;
    int fftIndex = blockIdx.y;
    int fa = fftIndex / (Nf1 * Nf2);
    int fb = (fftIndex / Nf2) % Nf1;
    int fc = fftIndex % Nf2;

    double * ms = mShifts + kPaddedDim * (fftIndex + Nf0*Nf1*Nf2 * mIndex);
    ms[0] = ma / double(Nr0 * Nf0) + fa / double(Nf0) + kgs0;
    ms[1] = mb / double(Nr1 * Nf1) + fb / double(Nf1) + kgs1;
    ms[2] = mc / double(Nr2 * Nf2) + fc / double(Nf2) + kgs2;
}

double * CuWhFourierTransformBase::createKshifts(const GeomVector3d &kFracGlobalShift) const
{
    double * d_kShifts = nullptr;
    unsigned size = kPaddedDim * Nr[0]*Nr[1]*Nr[2] * Nf[0]*Nf[1]*Nf[2];
    CUDA_CHECK ( cudaMalloc(reinterpret_cast<void **>(&d_kShifts), sizeof(double) * size) );
    int Nrep = Nr[0] * Nr[1] * Nr[2];
    int Nfft = Nf[0] * Nf[1] * Nf[2];
    int nThreadsPerBlock = 256;
    dim3 nBlocks( (Nrep+nThreadsPerBlock-1) / nThreadsPerBlock,
                  Nfft,
                  1
                 );
    populateKshifts<<<nBlocks, nThreadsPerBlock, 0, stream>>>
        (d_kShifts, Nr[0], Nr[1], Nr[2], Nf[0], Nf[1], Nf[2], kPaddedDim,
         kFracGlobalShift.at(0), kFracGlobalShift.at(1), kFracGlobalShift.at(2));
    return d_kShifts;
}

GeomVector3d CuWhFourierTransformBase::toFracK(const GeomVector3d &k_au) const
{
    GeomVector3d fracOffset( { k_au.dot(tb.latticeVectors[0]),
                               k_au.dot(tb.latticeVectors[1]),
                               k_au.dot(tb.latticeVectors[2]) });
    return fracOffset / (2 * std::numbers::pi);
}

/***************************************** HK *********************************************/

CuWhFourierTransformHk::CuWhFourierTransformHk(cudaStream_t stream, const TightBindingParameter_t &tb,
                                               const FourierTransformParameter_t &ftp, const std::array<unsigned, 3> &Noffset) :
    CuWhFourierTransformBase(stream, tb, ftp, Noffset, TB_OP_CNT)
{
    unsigned baseSize = sizeof(std::complex<double>) * Nf[0] * Nf[1] * Nf[2] * numWann * numWann;
    CUDA_CHECK( cudaMalloc(reinterpret_cast<void**>(&d_realSpaceMergedOps), baseSize * TB_OP_CNT  ) );
    CUDA_CHECK( cudaMemcpyAsync(d_realSpaceMergedOps + baseSize * TB_OP_POS(H0),
                                d_realSpaceGridOps + baseSize * TB_FULL_OP_POS(H0),
                                baseSize, cudaMemcpyDeviceToDevice, stream) );
}

CuWhFourierTransformHk::~CuWhFourierTransformHk()
{
    CUDA_CHECK ( cudaFree(d_realSpaceMergedOps) );
}

__device__ inline cuDoubleComplex cuCscale(double alpha, cuDoubleComplex v){
    return cuDoubleComplex( alpha * cuCreal(v), alpha * cuCimag(v) );
}

__global__ void mergeOps(const cuDoubleComplex * fullOps, cuDoubleComplex * mergedOps,
                    double E0, double E1, double E2, int baseSize)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if ( id >= baseSize )
        return;
    mergedOps[baseSize * TB_OP_POS(HE) + id] = cuCadd(fullOps[baseSize * TB_FULL_OP_POS(H0) + id],
                                               cuCadd(cuCscale(E0, fullOps[baseSize * (TB_FULL_OP_POS(D)+0) + id]),
                                               cuCadd(cuCscale(E1, fullOps[baseSize * (TB_FULL_OP_POS(D)+1) + id]),
                                               cuCscale(E2, fullOps[baseSize * (TB_FULL_OP_POS(D)+2) + id]))));
}

void CuWhFourierTransformHk::transform(cuDoubleComplex *data, const GeomVector3d &kOffset_au, const GeomVector3d &E)
{
    updateMeshShifts(toFracK(kOffset_au));
    unsigned baseSize = Nf[0] * Nf[1] * Nf[2] * numWann * numWann;
    unsigned nThreadsPerBlock = 32;
    mergeOps<<< (baseSize+nThreadsPerBlock-1)/nThreadsPerBlock, nThreadsPerBlock, 0, stream>>>
            (d_realSpaceGridOps, d_realSpaceMergedOps, E.at(0), E.at(1), E.at(2), baseSize);
    applyTransform(d_realSpaceMergedOps, data);
}


HkVec CuWhFourierTransformHk::unravelFFTmesh(const cuDoubleComplex * d_kOps) const // [ops][Nk][m][n]
{
    unsigned size = TB_OP_CNT * Nr[0]*Nr[1]*Nr[2] * Nf[0]*Nf[1]*Nf[2] * numWann * numWann;
    std::vector<std::complex<double>> data(size);
    CUDA_CHECK(cudaMemcpyAsync(data.data(), d_kOps, sizeof(std::complex<double>) * size,
                                 cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    std::array<unsigned, 3> Nk({Nf[0]*Nr[0], Nf[1]*Nr[1], Nf[2]*Nr[2]});
    HkVec res(numWann, Nk[0] * Nk[1] * Nk[2]);
    for(unsigned a=0; a<Nk[0]; a++)
        for(unsigned b=0; b<Nk[1]; b++)
            for(unsigned c=0; c<Nk[2]; c++){
                unsigned fa = a / Nr[0];
                unsigned fb = b / Nr[1];
                unsigned fc = c / Nr[2];
                unsigned outerId = (c % Nr[2]) + Nr[2] * ( (b%Nr[1]) + Nr[1] * (a%Nr[0]));
                unsigned oldId = fc + ftp.Nf[2] * (fb + ftp.Nf[1] * (fa + ftp.Nf[0] * outerId));
                unsigned newId = c + Nk[2] * (b + Nk[1] * a);
                for(unsigned matId=0; matId<numWann*numWann; matId++){
                    res.H0(newId)[matId] = data[ matId + numWann*numWann * (oldId + Nk[0]*Nk[1]*Nk[2] * TB_OP_POS(H0))];
                    res.HE(newId)[matId] = data[ matId + numWann*numWann * (oldId + Nk[0]*Nk[1]*Nk[2] * TB_OP_POS(HE))];
                }
            }
    return res;
}


/*************************************** HKfull *******************************************/

CuWhFourierTransformHkFull::CuWhFourierTransformHkFull(cudaStream_t stream, const TightBindingParameter_t &tb,
                                               const FourierTransformParameter_t &ftp, const std::array<unsigned, 3> &Noffset) :
    CuWhFourierTransformBase(stream, tb, ftp, Noffset, TB_FULL_OP_CNT)
{

}

CuWhFourierTransformHkFull::~CuWhFourierTransformHkFull()
{

}

void CuWhFourierTransformHkFull::transform(cuDoubleComplex *data, const GeomVector3d &kOffset_au)
{
    updateMeshShifts(toFracK(kOffset_au));
    applyTransform(d_realSpaceGridOps, data);
}


HkFullVec CuWhFourierTransformHkFull::unravelFFTmesh(const cuDoubleComplex * d_kOps) const // [ops][Nk][m][n]
{
    unsigned size = TB_FULL_OP_CNT * Nr[0]*Nr[1]*Nr[2] * Nf[0]*Nf[1]*Nf[2] * numWann * numWann;
    std::vector<std::complex<double>> data(size);
    CUDA_CHECK(cudaMemcpyAsync(data.data(), d_kOps, sizeof(std::complex<double>) * size,
                                 cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    std::array<unsigned, 3> Nk({Nf[0]*Nr[0], Nf[1]*Nr[1], Nf[2]*Nr[2]});
    HkFullVec res(numWann, Nk[0] * Nk[1] * Nk[2]);
    for(unsigned a=0; a<Nk[0]; a++)
        for(unsigned b=0; b<Nk[1]; b++)
            for(unsigned c=0; c<Nk[2]; c++){
                unsigned fa = a / Nr[0];
                unsigned fb = b / Nr[1];
                unsigned fc = c / Nr[2];
                unsigned outerId = (c % Nr[2]) + Nr[2] * ( (b%Nr[1]) + Nr[1] * (a%Nr[0]));
                unsigned oldId = fc + ftp.Nf[2] * (fb + ftp.Nf[1] * (fa + ftp.Nf[0] * outerId));
                unsigned newId = c + Nk[2] * (b + Nk[1] * a);
                for(unsigned matId=0; matId<numWann*numWann; matId++){
                    res.H0(newId)[matId] = data[ matId + numWann*numWann * (oldId + Nk[0]*Nk[1]*Nk[2] * TB_FULL_OP_POS(H0))];
                    for(unsigned d=0; d<3; d++){
                        res.D(newId, d)[matId] = data[ matId + numWann*numWann * (oldId + Nk[0]*Nk[1]*Nk[2] * (TB_FULL_OP_POS(D)+d))];
                        res.dH0dk(newId, d)[matId] = data[ matId + numWann*numWann * (oldId + Nk[0]*Nk[1]*Nk[2] * (TB_FULL_OP_POS(dH0dk)+d))];
                    }
                }
            }
    return res;
}
