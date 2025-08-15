#ifndef SBE_WH_FOURIER_TRANSFORM_H
#define SBE_WH_FOURIER_TRANSFORM_H


#include "util/w90util.h"
#include "util/GeomVector.hpp"

#include<vector>
#include<complex>

#include<numbers>

#include<fftw3.h>

#include "io/TightBindingParameter.h"
#include "io/FourierTransformParameter.h"


template<typename Tdata>
class fftw_allocator : public std::allocator<Tdata>
{
public:
    template <typename U>
    struct rebind { typedef fftw_allocator<U> other; };
    Tdata* allocate(size_t n) { return (Tdata*) fftw_malloc(sizeof(Tdata) * n); }
    void deallocate(Tdata* data, std::size_t size) { fftw_free(data); }
};


struct RealSpaceCell_t {
    CellIndex r;
    std::complex<double> H0;
    std::complex<double> D[3];
};

struct TbFullOperators_t {
    std::complex<double> H0;
    std::complex<double> D[3];
    std::complex<double> dH0dk[3];
};

#define TB_FULL_OP_CNT ( sizeof(TbFullOperators_t) / sizeof(std::complex<double>) )
#define TB_FULL_OP_POS(name) ( offsetof(TbFullOperators_t, name) / sizeof(std::complex<double>) )

struct TbOperators_t {
    std::complex<double> H0;
    std::complex<double> HE; // H + sum D_dir E_dir
};

#define TB_OP_CNT ( sizeof(TbOperators_t) / sizeof(std::complex<double>) )
#define TB_OP_POS(name) ( offsetof(TbOperators_t, name) / sizeof(std::complex<double>) )

class HkVec {
    const unsigned strideAlignment = 32;
    unsigned numWann;
    unsigned stride;
    unsigned numKpts;
    std::vector<std::complex<double>> d;
public:
    HkVec() { numWann = 0; numKpts = 0; }
    HkVec(unsigned numWann, unsigned numKpts) :
        numWann(numWann), stride( strideAlignment * (((numWann*numWann) + strideAlignment-1) / strideAlignment)),
        numKpts(numKpts), d( TB_OP_CNT * stride * numKpts, 0) {}

#define HK_VEC_INDEX(name, kid) ( stride * (TB_OP_CNT*kid + TB_OP_POS(name) ) )
#define HK_VEC_SG(name) \
    inline std::complex<double> * name(unsigned kid) { \
        return d.data() + HK_VEC_INDEX(name, kid); } \
    inline std::complex<double> get##name(unsigned kid, unsigned m, unsigned n) const { \
        return d[HK_VEC_INDEX(name, kid) + numWann * m + n]; }
    HK_VEC_SG(H0)
    HK_VEC_SG(HE)
#undef HK_VEC_SG
    inline void write(unsigned kid, unsigned m, unsigned n, TbOperators_t op){
#define HK_VEC_W(name) d[HK_VEC_INDEX(name, kid) + numWann * m + n] = op.name;
        HK_VEC_W(H0)
        HK_VEC_W(HE)
#undef HK_VEC_W
    }
#undef HK_VEC_INDEX
    inline std::complex<double> * data() { return d.data(); }
    inline unsigned getNumKpts() const { return numKpts; }
    inline unsigned getStride() const { return stride; }
};


class HkFullVec {
    const unsigned strideAlignment = 32;
    unsigned numWann;
    unsigned stride;
    unsigned numKpts;
    std::vector<std::complex<double>> d;
public:
    HkFullVec() { numWann = 0; numKpts = 0; }
    HkFullVec(unsigned numWann, unsigned numKpts) :
        numWann(numWann), stride( strideAlignment * (((numWann*numWann) + strideAlignment-1) / strideAlignment)),
        numKpts(numKpts), d(TB_FULL_OP_CNT * stride * numKpts, 0) {}
#define HK_FVEC_INDEX(name, kid) ( stride * (TB_FULL_OP_CNT*kid + TB_FULL_OP_POS(name) ) )
#define HK_FVEC_INDEXA(name, kid, dir) ( stride * (TB_FULL_OP_CNT*kid + TB_FULL_OP_POS(name) + dir ) )
#define HK_FVEC_SG(name) \
    inline std::complex<double> * name(unsigned kid) { \
        return d.data() + HK_FVEC_INDEX(name, kid); } \
    inline std::complex<double> get##name(unsigned kid, unsigned m, unsigned n) const { \
        return d[HK_FVEC_INDEX(name, kid) + numWann * m + n]; }
#define HK_FVEC_SGA(name) \
    inline std::complex<double> * name(unsigned kid, unsigned dir) { \
        return d.data() + HK_FVEC_INDEXA(name, kid, dir); } \
    inline std::complex<double> get##name(unsigned kid, unsigned dir, unsigned m, unsigned n) const { \
        return d[HK_FVEC_INDEXA(name, kid, dir) + numWann * m + n]; }
    HK_FVEC_SG(H0)
    HK_FVEC_SGA(D)
    HK_FVEC_SGA(dH0dk)
#undef HK_VEC_SG
#undef HK_VEC_SGA

    inline void write(unsigned kid, unsigned m, unsigned n, TbFullOperators_t op){
#define HK_FVEC_W(name) d[HK_FVEC_INDEX(name, kid) + numWann * m + n] = op.name;
#define HK_FVEC_WA(name, dir) d[HK_FVEC_INDEXA(name, kid, dir) + numWann * m + n] = op.name[dir];
        HK_FVEC_W(H0)
        HK_FVEC_WA(D, 0)
        HK_FVEC_WA(D, 1)
        HK_FVEC_WA(D, 2)
        HK_FVEC_WA(dH0dk, 0)
        HK_FVEC_WA(dH0dk, 1)
        HK_FVEC_WA(dH0dk, 2)
#undef HK_FVEC_W
#undef HK_FVEC_WA
    }
    inline std::complex<double> * data() { return d.data(); }
    inline unsigned getNumKpts() const { return numKpts; }
    inline unsigned getStride() const { return stride; }
};

class WhFourierTransform{
    const TightBindingParameter_t &tb;
    const FourierTransformParameter_t &ftp;
    unsigned numWann;
    unsigned dim;
    unsigned numCells;
    GeomVector3d kBase[3];
    std::array<unsigned, 3> Nrepeat;
    std::array<unsigned, 3> Nk;
    std::vector<GeomVector3d> shifts;
    std::vector<std::complex<double>, fftw_allocator<std::complex<double> > > fftData;
    std::vector<TbFullOperators_t> realDataGrid;
    std::vector<RealSpaceCell_t> realData;
    inline GeomVector3d toFracK(const GeomVector3d &k_au) const;
    inline unsigned rdIndex(unsigned m, unsigned n, unsigned cellIndex) const { return numCells * (m*numWann+n) + cellIndex; }
    inline unsigned gRelativeIndex(unsigned m, unsigned n, int a, int b, int c) const
                         { return c - tb.minCellIndices.at(2) + ftp.Nf[2] *
                                ( b - tb.minCellIndices.at(1) + ftp.Nf[1] *
                                ( a - tb.minCellIndices.at(0) + ftp.Nf[0] *
                                ( n + tb.numWann * m))); }

    inline unsigned gAbsoluteIndex(unsigned m, unsigned n, unsigned a, unsigned b, unsigned c) const
                         { return c + ftp.Nf[2] * (b + ftp.Nf[1] * (a + ftp.Nf[0] * ( n + numWann * m))); }

    inline unsigned fftIndex(unsigned outerId, unsigned m, unsigned n, unsigned opNr, unsigned a, unsigned b, unsigned c) const
             { return c + ftp.Nf[2] * (b + ftp.Nf[1] * (a + ftp.Nf[0] *
                         ( opNr + TB_OP_CNT * (n + numWann * (m + numWann * outerId))))); }
    inline unsigned fftIndex(unsigned m, unsigned n, unsigned opNr, unsigned a, unsigned b, unsigned c) const
             { return c + ftp.Nf[2] * (b + ftp.Nf[1] * (a + ftp.Nf[0] *
                         ( opNr + TB_OP_CNT * (n + numWann * m)))); }

    inline unsigned fftFullIndex(unsigned outerId, unsigned m, unsigned n, unsigned opNr, unsigned a, unsigned b, unsigned c) const
             { return c + ftp.Nf[2] * (b + ftp.Nf[1] * (a + ftp.Nf[0] *
                         ( opNr + TB_FULL_OP_CNT * (n + numWann * (m + numWann * outerId))))); }
    inline unsigned fftFullIndex(unsigned m, unsigned n, unsigned opNr, unsigned a, unsigned b, unsigned c) const
             { return c + ftp.Nf[2] * (b + ftp.Nf[1] * (a + ftp.Nf[0] *
                         ( opNr + TB_FULL_OP_CNT * (n + numWann *m)))); }

    TbFullOperators_t transformFracPoint(const GeomVector3d &kFrac, unsigned m, unsigned n) const;
    TbOperators_t transformFracPoint(const GeomVector3d &kFrac, unsigned m, unsigned n, const GeomVector3d &E) const;
    // for FFT
    fftw_plan fftPlan;
    fftw_plan fftFullPlan;
    unsigned offsetCount;
public:
    WhFourierTransform(const TightBindingParameter_t &tb, const FourierTransformParameter_t &ftp, std::array<unsigned, 3> Nk, unsigned dim);
    ~WhFourierTransform();
    std::vector<GeomVector3d> getKmesh(const GeomVector3d &kOffset_au) const;

    const std::vector<GeomVector3d> & getFineMeshFracShifts(const GeomVector3d &kOffset_au);
    const std::vector<GeomVector3d> & getFineMeshFracShifts(const GeomVector3d &kOffset_au, unsigned startIndex, unsigned count);
    GeomVector3d getFFTfracShift(unsigned kpt, const std::vector<GeomVector3d> &meshShifts, const GeomVector3d &externalShift_au);
    HkVec transformMeshDirect(const GeomVector3d &kOffset_au, const GeomVector3d &E) const;
    HkFullVec transformMeshDirect(const GeomVector3d &kOffset_au) const;
    bool planFFT(unsigned newOffsetCount) { return planFFT(newOffsetCount, ftp.fftwPlanningFlag, ftp.planningTimeLimit); }
    bool planFFT(unsigned newOffsetCount, unsigned fftwPlaningFlag, double timelimit = -1);
    HkVec createHkVec() const  { return HkVec(numWann, offsetCount * ftp.Nf[0] * ftp.Nf[1] * ftp.Nf[2] ); }
    HkFullVec createHkFullVec() const  { return HkFullVec(numWann, offsetCount * ftp.Nf[0] * ftp.Nf[1] * ftp.Nf[2] ); }

    void transformMeshFFT(HkVec & fftOps, const GeomVector3d &kFracOffset, const GeomVector3d &E);
    void transformMeshFFT(HkVec & fftOps, const std::vector<GeomVector3d> &kFracOffsets, const GeomVector3d &E);
    void transformMeshFFT(HkFullVec & fftFullOps, const GeomVector3d &kFracOffset);
    void transformMeshFFT(HkFullVec & fftFullOps, const std::vector<GeomVector3d> &kFracOffsets);
    HkVec unravelFFTmesh(const HkVec & fftOps) const;
    HkFullVec unravelFFTfullMesh(const HkFullVec & fftFullOps) const;
    void finalizeFFT();
};

#endif
