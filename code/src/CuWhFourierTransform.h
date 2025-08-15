#ifndef SBE_WH_CU_FOURIER_TRANSFORM_H
#define SBE_WH_CU_FOURIER_TRANSFORM_H

#include<array>
#include<vector>

#include <cuda_runtime.h>
#include <cufft.h>
#include "cuUtils.h"

#include "WhFourierTransform.h"

class CuWhFourierTransformBase {
protected:
    const unsigned kPaddedDim = 4;
    cudaStream_t stream;
    const TightBindingParameter_t &tb;
    const FourierTransformParameter_t &ftp;
    const std::array<unsigned, 3> &Nr;
    const std::array<unsigned, 3> &Nf;
    unsigned numOps;
    const CellIndex &minCellIndices;
    unsigned numWann;

    double * d_meshShifts;
    cuDoubleComplex * d_realSpaceGridOps; // [op][m,n][R]
    cuDoubleComplex * d_fftData;
    cufftHandle plan;

    void updateMeshShifts(const GeomVector3d &kFracOffset);
    // out: [op][Nr][Nf_k][m][n] ; in: [op][m][n][[R]
    void applyTransform(const cuDoubleComplex * rOps, cuDoubleComplex *kOps);
public:
    CuWhFourierTransformBase(cudaStream_t stream, const TightBindingParameter_t &tb, const FourierTransformParameter_t &ftp,
                             const std::array<unsigned, 3> &Nr, unsigned numOps);
    ~CuWhFourierTransformBase();
    GeomVector3d toFracK(const GeomVector3d &k_au) const;
    std::vector<GeomVector3d> getKmesh(const GeomVector3d &kOffset_au) const;
    double * createKshifts(const GeomVector3d &kFracGlobalShift) const; // on GPU
    cuDoubleComplex * createDataVec() const; // on GPU
    unsigned getOpCount() const { return numOps; }
    unsigned getKpaddedDim() const { return kPaddedDim; }
};

class CuWhFourierTransformHk : public CuWhFourierTransformBase {
    cuDoubleComplex * d_realSpaceMergedOps;
public:
    CuWhFourierTransformHk(cudaStream_t stream, const TightBindingParameter_t &tb, const FourierTransformParameter_t &ftp,
                           const std::array<unsigned, 3> &Noffset);
    ~CuWhFourierTransformHk();
    void transform(cuDoubleComplex *data, const GeomVector3d &kOffset_au, const GeomVector3d &E);

    HkVec unravelFFTmesh(const cuDoubleComplex * kOps) const; // [ops][Nk][m][n]
};

class CuWhFourierTransformHkFull : public CuWhFourierTransformBase {
public:
    CuWhFourierTransformHkFull(cudaStream_t stream, const TightBindingParameter_t &tb, const FourierTransformParameter_t &ftp,
                           const std::array<unsigned, 3> &Noffset);
    ~CuWhFourierTransformHkFull();
    void transform(cuDoubleComplex *data, const GeomVector3d &kOffset_au);
    HkFullVec unravelFFTmesh(const cuDoubleComplex * kOps) const; // [ops][Nk][m][n]
};


#endif
