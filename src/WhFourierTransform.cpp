#include "WhFourierTransform.h"

WhFourierTransform::WhFourierTransform(const TightBindingParameter_t &tb, const FourierTransformParameter_t &ftp, std::array<unsigned, 3> Nk, unsigned dim)
    : tb(tb), ftp(ftp), Nk(Nk), dim(dim)
{
    numCells = tb.realSpaceOps.size();
    numWann = tb.numWann;
    realData.resize(numCells * numWann * numWann);
    realDataGrid.resize(numWann * numWann * ftp.Nf[0] * ftp.Nf[1] * ftp.Nf[2], {0, 0, 0, 0});
    const auto &Rl = tb.latticeVectors;
    unsigned cid = 0;
    for(const auto [ci, Hmat] : tb.realSpaceOps){
        for(unsigned n=0; n<numWann; n++)
            for(unsigned m=0; m<numWann; m++){
                unsigned rid = rdIndex(m, n, cid);
                realData[rid].r = ci;
                GeomVector3d pos = (double)ci.at(0) * Rl[0] + (double)ci.at(1) * Rl[1] + (double)ci.at(2) * Rl[2];
                realData[rid].H0 = Hmat.H(m, n);
                for(unsigned dir=0; dir<3; dir++)
                    realData[rid].D[dir] = Hmat.D[dir](m, n);
                unsigned gid = gRelativeIndex(m, n, ci.at(0), ci.at(1), ci.at(2));
                realDataGrid[gid].H0 = Hmat.H(m, n);
                for(unsigned dir=0; dir<3; dir++){
                    realDataGrid[gid].D[dir] = Hmat.D[dir](m, n);
                    realDataGrid[gid].dH0dk[dir] = std::complex<double>{0, 1} * pos[dir]  * Hmat.H(m, n);
                }
            }
        cid++;
    }
    for(unsigned dir=0; dir<3; dir++)
        Nrepeat[dir] = Nk[dir]/ftp.Nf[dir];
    shifts.resize(Nrepeat[0] * Nrepeat[1] * Nrepeat[2]);
    offsetCount = 0;
    double preFact = 1  / ( Rl[0].dot(Rl[1].cross(Rl[2])));
    kBase[0] = preFact * Rl[1].cross(Rl[2]);
    kBase[1] = preFact * Rl[2].cross(Rl[0]);
    kBase[2] = preFact * Rl[0].cross(Rl[1]);
}

WhFourierTransform::~WhFourierTransform()
{
    if (offsetCount)
        finalizeFFT();
}

GeomVector3d WhFourierTransform::toFracK(const GeomVector3d &k_au) const
{
    GeomVector3d fracOffset( { k_au.dot(tb.latticeVectors[0]),
                               k_au.dot(tb.latticeVectors[1]),
                               k_au.dot(tb.latticeVectors[2]) });
    return fracOffset / (2 * std::numbers::pi);
}

TbOperators_t WhFourierTransform::transformFracPoint(const GeomVector3d &kFrac, unsigned m, unsigned n,
                                                     const GeomVector3d &E) const
{
    TbOperators_t Hk = {};
    for(unsigned i=0; i<numCells; i++){
        const RealSpaceCell_t & rc = realData[rdIndex(m, n, i)];
        double phase = 2 * std::numbers::pi * (kFrac.at(0) * rc.r.at(0) + kFrac.at(1) * rc.r.at(1) + kFrac.at(2) * rc.r.at(2));
        std::complex<double> eikr(cos(phase), sin(phase));
        Hk.H0 += rc.H0 * eikr;
        std::complex<double> rHE = rc.H0;
        for(unsigned i=0; i<dim; i++)
            rHE += rc.D[i] * E.at(i);
        Hk.HE += rHE * eikr;
    }
    return Hk;
}

TbFullOperators_t WhFourierTransform::transformFracPoint(const GeomVector3d &kFrac, unsigned m, unsigned n) const
{
    TbFullOperators_t Hk = {};
    for(unsigned ci=0; ci<numCells; ci++){
        const RealSpaceCell_t & rc = realData[rdIndex(m, n, ci)];
        double phase = 2 * std::numbers::pi * (kFrac.at(0) * rc.r.at(0) + kFrac.at(1) * rc.r.at(1) + kFrac.at(2) * rc.r.at(2));
        std::complex<double> eikr(cos(phase), sin(phase));
        Hk.H0 += rc.H0 * eikr;
        for(unsigned i=0; i<dim; i++)
            Hk.D[i] += rc.D[i] * eikr;
    }
    return Hk;
}


std::vector<GeomVector3d> WhFourierTransform::getKmesh(const GeomVector3d &kOffset_au) const
{
    GeomVector3d kFracOffset = toFracK(kOffset_au);
    std::vector<GeomVector3d> kMesh;
    kMesh.reserve(Nk[0] * Nk[1] * Nk[2]);
    for(unsigned a=0; a<Nk[0]; a++)
        for(unsigned b=0; b<Nk[1]; b++)
            for(unsigned c=0; c<Nk[2]; c++)
                kMesh.push_back( a /(double)Nk[0] * kBase[0] + b/(double)Nk[1] * kBase[1] + c/(double)Nk[2] * kBase[2]
                                 + kFracOffset );
    return kMesh;
}

HkVec WhFourierTransform::transformMeshDirect(const GeomVector3d &kOffset_au, const GeomVector3d &E) const
{
    HkVec res(numWann, Nk[0] * Nk[1] * Nk[2]);
    GeomVector3d kFracOffset = toFracK(kOffset_au);
    #pragma omp parallel for
    for(unsigned id=0; id<Nk[0]*Nk[1]*Nk[2]; id++){
        for(unsigned m=0; m<numWann; m++)
            for(unsigned n=0; n<numWann; n++){
                unsigned a = id / (Nk[1]*Nk[2]);
                unsigned b = (id / Nk[2]) % Nk[1];
                unsigned c = id % Nk[2];
                GeomVector3d k( {a /(double)Nk[0], b/(double)Nk[1],  c/(double)Nk[2]});
                k += kFracOffset;
                res.write(id, m, n, transformFracPoint(k, m, n, E));
            }
    }
    return res;
}

HkFullVec WhFourierTransform::transformMeshDirect(const GeomVector3d &kOffset_au) const
{
    HkFullVec res(numWann, Nk[0] * Nk[1] * Nk[2]);
    GeomVector3d kFracOffset = toFracK(kOffset_au);
    #pragma omp parallel for
    for(unsigned id=0; id<Nk[0]*Nk[1]*Nk[2]; id++){
        for(unsigned m=0; m<numWann; m++)
            for(unsigned n=0; n<numWann; n++){
                unsigned a = id / (Nk[1]*Nk[2]);
                unsigned b = (id / Nk[2]) % Nk[1];
                unsigned c = id % Nk[2];
                GeomVector3d k( {a /(double)Nk[0], b/(double)Nk[1],  c/(double)Nk[2]});
                k += kFracOffset;
                res.write(id, m, n, transformFracPoint(k, m, n));
            }
    }
    return res;
}

const std::vector<GeomVector3d> & WhFourierTransform::getFineMeshFracShifts(const GeomVector3d &kOffset_au, unsigned startIndex, unsigned count)
{
    assert ( count == shifts.size() );
    GeomVector3d kFracOffset = toFracK(kOffset_au);
    unsigned cc = 0;
    for(unsigned a=startIndex/(Nrepeat[1]*Nrepeat[2]); a<Nrepeat[0] && cc < count; a++)
        for(unsigned b=(startIndex/Nrepeat[2])%Nrepeat[1]; b<Nrepeat[1] && cc < count; b++)
            for(unsigned c=startIndex%Nrepeat[2]; c<Nrepeat[2] && cc < count; c++){
                GeomVector3d k({ a/(double)Nk[0], b/(double)Nk[1], c/(double)Nk[2] });
                shifts[cc++] = k + kFracOffset;
            }
    return shifts;
}

const std::vector<GeomVector3d> & WhFourierTransform::getFineMeshFracShifts(const GeomVector3d &kOffset_au)
{
    GeomVector3d kFracOffset = toFracK(kOffset_au);
    for(unsigned a=0; a<Nrepeat[0]; a++)
        for(unsigned b=0; b<Nrepeat[1]; b++)
            for(unsigned c=0; c<Nrepeat[2]; c++){
                GeomVector3d k({ a/(double)Nk[0], b/(double)Nk[1], c/(double)Nk[2] });
                shifts[a*Nrepeat[1]*Nrepeat[2] + b*Nrepeat[2] + c] = k + kFracOffset;
            }
    return shifts;
}

GeomVector3d WhFourierTransform::getFFTfracShift(unsigned kpt, const std::vector<GeomVector3d> &meshShifts, const GeomVector3d &externalShift_au)
{
    unsigned Nf[3] = { Nk[0] / Nrepeat[0], Nk[1] / Nrepeat[1], Nk[2] / Nrepeat[2]};
    unsigned meshId = kpt / (Nf[0] * Nf[1] * Nf[2]);
    unsigned fftId = kpt % (Nf[0] * Nf[1] * Nf[2]);
    GeomVector3d kShift = meshShifts[meshId] + toFracK(externalShift_au);
    unsigned a = fftId / (Nf[1]*Nf[2]);
    unsigned b = (fftId / Nf[2]) % Nf[1];
    unsigned c = fftId % Nf[2];
    GeomVector3d meshShift( {a /(double)Nf[0], b/(double)Nf[1],  c/(double)Nf[2]});
    return kShift + meshShift;
}

bool WhFourierTransform::planFFT(unsigned newOffsetCount, unsigned fftwPlanningFlag, double timeLimit)
{
    if ( offsetCount == 0 && offsetCount == newOffsetCount )
        return false;
    if ( offsetCount == newOffsetCount )
        return true;
    if ( offsetCount != 0)
        finalizeFFT();
    offsetCount = newOffsetCount;
    int dist = ftp.Nf[0] * ftp.Nf[1] * ftp.Nf[2];
    assert(TB_FULL_OP_CNT >= TB_OP_CNT);
    fftData.resize(TB_FULL_OP_CNT * offsetCount * dist * numWann * numWann, 0);
    int dims[3] = { (int)ftp.Nf[0], (int)ftp.Nf[1], (int)ftp.Nf[2] };
    fftw_set_timelimit(timeLimit);
    // transform inner most indices, i.e.  k_x, k_y, k_z
    // remaining indices ordered (from inner) to (outer) {H0|He}, m, n, offsetCount
    fftPlan = fftw_plan_many_dft(3, dims, (int)offsetCount * numWann * numWann * TB_OP_CNT,
                                 reinterpret_cast<fftw_complex*>(&fftData[0]), NULL,
                                 1, dist,
                                 reinterpret_cast<fftw_complex*>(&fftData[0]), NULL,
                                 1, dist,
                                 FFTW_BACKWARD, fftwPlanningFlag);
    fftFullPlan = fftw_plan_many_dft(3, dims, (int)offsetCount * numWann * numWann * TB_FULL_OP_CNT,
                                 reinterpret_cast<fftw_complex*>(&fftData[0]), NULL,
                                 1, dist,
                                 reinterpret_cast<fftw_complex*>(&fftData[0]), NULL,
                                 1, dist,
                                 FFTW_BACKWARD, fftwPlanningFlag);
    shifts.resize(newOffsetCount);
    return true;
}

void WhFourierTransform::transformMeshFFT(HkVec &fftOps,
                                          const GeomVector3d &kFracOffset,
                                          const GeomVector3d &E)
{
    assert(offsetCount == 1);
    // shift phases -- kOffset
    #pragma omp parallel for
    for(unsigned id=0; id<numWann * numWann; id++){
        unsigned m = id / numWann;
        unsigned n = id % numWann;
        for(int a=0; a<ftp.Nf[0]; a++)
            for(int b=0; b<ftp.Nf[1]; b++)
                for(int c=0; c<ftp.Nf[2]; c++){
                    double phase = 2 * std::numbers::pi * ( a * kFracOffset.at(0) +
                                                            b * kFracOffset.at(1) +
                                                            c * kFracOffset.at(2));
                    unsigned gi = gAbsoluteIndex(m, n, a, b, c);
                    std::complex<double> eikr(cos(phase), sin(phase));
                    fftData[fftIndex(m, n, TB_OP_POS(H0), a, b, c)] = realDataGrid[gi].H0 * eikr;
                    std::complex<double> re = realDataGrid[gi].H0;
                    for(unsigned dir=0; dir<dim; dir++)
                        re += realDataGrid[gi].D[dir] * E.at(dir);
                    fftData[fftIndex(m, n, TB_OP_POS(HE), a, b, c)] = re * eikr;
                }
    }
    // apply FFTs
    fftw_execute(fftPlan);
    // shift phases -- rOffset
    double rPhaseFrac[3] = { tb.minCellIndices.at(0) / (double)ftp.Nf[0], tb.minCellIndices.at(1) / (double)ftp.Nf[1],
                       tb.minCellIndices.at(2) / (double)ftp.Nf[2] };
    #pragma omp parallel for
    for(unsigned id=0; id<numWann * numWann; id++){
        double offFracPhase = tb.minCellIndices.at(0) * kFracOffset.at(0) +
                              tb.minCellIndices.at(1) * kFracOffset.at(1) +
                              tb.minCellIndices.at(2) * kFracOffset.at(2);
        unsigned m = id / numWann;
        unsigned n = id % numWann;
        for(int a=0; a<ftp.Nf[0]; a++)
            for(int b=0; b<ftp.Nf[1]; b++)
                for(int c=0; c<ftp.Nf[2]; c++){
                    double phase = 2 * std::numbers::pi * ( offFracPhase +
                            a * rPhaseFrac[0] + b * rPhaseFrac[1] + c * rPhaseFrac[2]);
                    unsigned gi = gAbsoluteIndex(m, n, a, b, c);
                    std::complex<double> eikr(cos(phase), sin(phase));
                    unsigned oldId = c + ftp.Nf[2] * (b + ftp.Nf[1] * a);
                    fftOps.H0(oldId)[id] = fftData[fftIndex(m, n, TB_OP_POS(H0), a, b, c)] * eikr;
                    fftOps.HE(oldId)[id] = fftData[fftIndex(m, n, TB_OP_POS(HE), a, b, c)] * eikr;
                }
    }
}

void WhFourierTransform::transformMeshFFT(HkVec &fftOps,
                                          const std::vector<GeomVector3d> &kFracOffsets,
                                          const GeomVector3d &E)
{
    assert(offsetCount == kFracOffsets.size());
    if ( kFracOffsets.size() == 1 ){
        transformMeshFFT(fftOps, kFracOffsets[0], E);
        return;
    }
    // shift phases -- kOffset
    #pragma omp parallel for
    for(unsigned id=0; id<offsetCount * numWann * numWann; id++){
        unsigned oid = id / (numWann * numWann);
        unsigned m = (id %(numWann*numWann)) / numWann;
        unsigned n = id % numWann;
        for(int a=0; a<ftp.Nf[0]; a++)
            for(int b=0; b<ftp.Nf[1]; b++)
                for(int c=0; c<ftp.Nf[2]; c++){
                    double phase = 2 * std::numbers::pi * ( a * kFracOffsets[oid].at(0) +
                                                            b * kFracOffsets[oid].at(1) +
                                                            c * kFracOffsets[oid].at(2));
                    unsigned gi = gAbsoluteIndex(m, n, a, b, c);
                    std::complex<double> eikr(cos(phase), sin(phase));
                    fftData[fftIndex(oid, m, n, TB_OP_POS(H0), a, b, c)] = realDataGrid[gi].H0 * eikr;
                    std::complex<double> re = realDataGrid[gi].H0;
                    for(unsigned dir=0; dir<dim; dir++)
                        re += realDataGrid[gi].D[dir] * E.at(dir);
                    fftData[fftIndex(oid, m, n, TB_OP_POS(HE), a, b, c)] = re * eikr;
                }
    }
    // apply FFTs
    fftw_execute(fftPlan);
    // shift phases -- rOffset
    double rPhaseFrac[3] = { tb.minCellIndices.at(0) / (double)ftp.Nf[0], tb.minCellIndices.at(1) / (double)ftp.Nf[1],
                       tb.minCellIndices.at(2) / (double)ftp.Nf[2] };
    #pragma omp parallel for
    for(unsigned id=0; id<offsetCount * numWann * numWann; id++){
        unsigned oid = id / (numWann * numWann);
        double offFracPhase = tb.minCellIndices.at(0) * kFracOffsets[oid].at(0) +
                              tb.minCellIndices.at(1) * kFracOffsets[oid].at(1) +
                              tb.minCellIndices.at(2) * kFracOffsets[oid].at(2);
        unsigned matIndex = id%(numWann*numWann);
        unsigned m = matIndex / numWann;
        unsigned n = id % numWann;
        for(int a=0; a<ftp.Nf[0]; a++)
            for(int b=0; b<ftp.Nf[1]; b++)
                for(int c=0; c<ftp.Nf[2]; c++){
                    double phase = 2 * std::numbers::pi * ( offFracPhase +
                            a * rPhaseFrac[0] + b * rPhaseFrac[1] + c * rPhaseFrac[2]);
                    unsigned gi = gAbsoluteIndex(m, n, a, b, c);
                    std::complex<double> eikr(cos(phase), sin(phase));
                    unsigned oldId = c + ftp.Nf[2] * (b + ftp.Nf[1] * (a + ftp.Nf[0] * oid));
                    fftOps.H0(oldId)[matIndex] = fftData[fftIndex(oid, m, n, TB_OP_POS(H0), a, b, c)] * eikr;
                    fftOps.HE(oldId)[matIndex] = fftData[fftIndex(oid, m, n, TB_OP_POS(HE), a, b, c)] * eikr;
                }
    }
}

void WhFourierTransform::transformMeshFFT(HkFullVec &fftFullOps, const GeomVector3d &kFracOffset)
{
    assert(offsetCount == 1);
    // shift phases -- kOffset
    #pragma omp parallel for
    for(unsigned id=0; id<numWann * numWann; id++){
        unsigned m = id / numWann;
        unsigned n = id % numWann;
        for(int a=0; a<ftp.Nf[0]; a++)
            for(int b=0; b<ftp.Nf[1]; b++)
                for(int c=0; c<ftp.Nf[2]; c++){
                    double phase = 2 * std::numbers::pi * ( a * kFracOffset.at(0) +
                                                            b * kFracOffset.at(1) +
                                                            c * kFracOffset.at(2));
                    unsigned gi = gAbsoluteIndex(m, n, a, b, c);
                    std::complex<double> eikr(cos(phase), sin(phase));
                    fftData[fftFullIndex(m, n, TB_FULL_OP_POS(H0), a, b, c)] = realDataGrid[gi].H0 * eikr;
                    std::complex<double> re = realDataGrid[gi].H0;
                    for(unsigned dir=0; dir<dim; dir++){
                        fftData[fftFullIndex(m, n, TB_FULL_OP_POS(D)+dir, a, b, c)] = realDataGrid[gi].D[dir] * eikr;
                        fftData[fftFullIndex(m, n, TB_FULL_OP_POS(dH0dk)+dir, a, b, c)] = realDataGrid[gi].dH0dk[dir] * eikr;
                    }
                }
    }
    // apply FFTs
    fftw_execute(fftFullPlan);
    // shift phases -- rOffset
    double rPhaseFrac[3] = { tb.minCellIndices.at(0) / (double)ftp.Nf[0], tb.minCellIndices.at(1) / (double)ftp.Nf[1],
                       tb.minCellIndices.at(2) / (double)ftp.Nf[2] };
    double offFracPhase = tb.minCellIndices.at(0) * kFracOffset.at(0) +
                          tb.minCellIndices.at(1) * kFracOffset.at(1) +
                          tb.minCellIndices.at(2) * kFracOffset.at(2);
    #pragma omp parallel for
    for(unsigned id=0; id<numWann * numWann; id++){
        unsigned m = id / numWann;
        unsigned n = id % numWann;
        for(int a=0; a<ftp.Nf[0]; a++)
            for(int b=0; b<ftp.Nf[1]; b++)
                for(int c=0; c<ftp.Nf[2]; c++){
                    double phase = 2 * std::numbers::pi * ( offFracPhase +
                            a * rPhaseFrac[0] + b * rPhaseFrac[1] + c * rPhaseFrac[2]);
                    unsigned gi = gAbsoluteIndex(m, n, a, b, c);
                    std::complex<double> eikr(cos(phase), sin(phase));
                    unsigned oldId = c + ftp.Nf[2] * (b + ftp.Nf[1] * a);
                    fftFullOps.H0(oldId)[id] = fftData[fftFullIndex(m, n, TB_FULL_OP_POS(H0), a, b, c)] * eikr;
                    for(unsigned dir=0; dir<3; dir++){
                        fftFullOps.D(oldId, dir)[id] = fftData[fftFullIndex(m, n, TB_FULL_OP_POS(D)+dir, a, b, c)] * eikr;
                        fftFullOps.dH0dk(oldId, dir)[id] = fftData[fftFullIndex(m, n, TB_FULL_OP_POS(dH0dk)+dir, a, b, c)] * eikr;
                    }
                }
    }
}

void WhFourierTransform::transformMeshFFT(HkFullVec &fftFullOps, const std::vector<GeomVector3d> &kFracOffsets)
{
    assert(offsetCount == kFracOffsets.size());
    if ( kFracOffsets.size() == 1){
        transformMeshFFT(fftFullOps, kFracOffsets[0]);
        return;
    }
    // shift phases -- kOffset
    #pragma omp parallel for
    for(unsigned id=0; id<offsetCount * numWann * numWann; id++){
        unsigned oid = id / (numWann * numWann);
        unsigned m = (id %(numWann*numWann)) / numWann;
        unsigned n = id % numWann;
        for(int a=0; a<ftp.Nf[0]; a++)
            for(int b=0; b<ftp.Nf[1]; b++)
                for(int c=0; c<ftp.Nf[2]; c++){
                    double phase = 2 * std::numbers::pi * ( a * kFracOffsets[oid].at(0) +
                                                            b * kFracOffsets[oid].at(1) +
                                                            c * kFracOffsets[oid].at(2));
                    unsigned gi = gAbsoluteIndex(m, n, a, b, c);
                    std::complex<double> eikr(cos(phase), sin(phase));
                    fftData[fftFullIndex(oid, m, n, TB_FULL_OP_POS(H0), a, b, c)] = realDataGrid[gi].H0 * eikr;
                    std::complex<double> re = realDataGrid[gi].H0;
                    for(unsigned dir=0; dir<dim; dir++){
                        fftData[fftFullIndex(oid, m, n, TB_FULL_OP_POS(D)+dir, a, b, c)] = realDataGrid[gi].D[dir] * eikr;
                        fftData[fftFullIndex(oid, m, n, TB_FULL_OP_POS(dH0dk)+dir, a, b, c)] = realDataGrid[gi].dH0dk[dir] * eikr;
                    }
                }
    }
    // apply FFTs
    fftw_execute(fftFullPlan);
    // shift phases -- rOffset
    double rPhaseFrac[3] = { tb.minCellIndices.at(0) / (double)ftp.Nf[0], tb.minCellIndices.at(1) / (double)ftp.Nf[1],
                       tb.minCellIndices.at(2) / (double)ftp.Nf[2] };
    #pragma omp parallel for
    for(unsigned id=0; id<offsetCount * numWann * numWann; id++){
        unsigned oid = id / (numWann * numWann);
        double offFracPhase = tb.minCellIndices.at(0) * kFracOffsets[oid].at(0) +
                              tb.minCellIndices.at(1) * kFracOffsets[oid].at(1) +
                              tb.minCellIndices.at(2) * kFracOffsets[oid].at(2);
        unsigned matIndex = id%(numWann*numWann);
        unsigned m = matIndex / numWann;
        unsigned n = id % numWann;
        for(int a=0; a<ftp.Nf[0]; a++)
            for(int b=0; b<ftp.Nf[1]; b++)
                for(int c=0; c<ftp.Nf[2]; c++){
                    double phase = 2 * std::numbers::pi * ( offFracPhase +
                            a * rPhaseFrac[0] + b * rPhaseFrac[1] + c * rPhaseFrac[2]);
                    unsigned gi = gAbsoluteIndex(m, n, a, b, c);
                    std::complex<double> eikr(cos(phase), sin(phase));
                    unsigned oldId = c + ftp.Nf[2] * (b + ftp.Nf[1] * (a + ftp.Nf[0] * oid));
                    fftFullOps.H0(oldId)[matIndex] = fftData[fftFullIndex(oid, m, n, TB_FULL_OP_POS(H0), a, b, c)] * eikr;
                    for(unsigned dir=0; dir<3; dir++){
                        fftFullOps.D(oldId, dir)[matIndex] = fftData[fftFullIndex(oid, m, n, TB_FULL_OP_POS(D)+dir, a, b, c)] * eikr;
                        fftFullOps.dH0dk(oldId, dir)[matIndex] = fftData[fftFullIndex(oid, m, n, TB_FULL_OP_POS(dH0dk)+dir, a, b, c)] * eikr;
                    }
                }
    }
}

void WhFourierTransform::finalizeFFT()
{
    fftw_destroy_plan(fftPlan);
    fftw_destroy_plan(fftFullPlan);
}

HkVec WhFourierTransform::unravelFFTmesh(const HkVec &fftOps) const {
    assert( fftOps.getNumKpts() == ftp.Nf[0] * ftp.Nf[1] * ftp.Nf[2] * Nrepeat[0] * Nrepeat[1] * Nrepeat[2]);
    HkVec res(numWann, Nk[0] * Nk[1] * Nk[2]);
    for(unsigned a=0; a<Nk[0]; a++)
        for(unsigned b=0; b<Nk[1]; b++)
            for(unsigned c=0; c<Nk[2]; c++){
                unsigned fa = a / Nrepeat[0];
                unsigned fb = b / Nrepeat[1];
                unsigned fc = c / Nrepeat[2];
                unsigned outerId = (c % Nrepeat[2]) + Nrepeat[2] * ( (b%Nrepeat[1]) + Nrepeat[1] * (a%Nrepeat[0]));
                unsigned oldId = fc + ftp.Nf[2] * (fb + ftp.Nf[1] * (fa + ftp.Nf[0] * outerId));
                unsigned newId = c + Nk[2] * (b + Nk[1] * a);
                for(unsigned matId=0; matId<numWann*numWann; matId++){
                    res.H0(newId)[matId] = fftOps.getH0(oldId, matId/numWann, matId % numWann);
                    res.HE(newId)[matId] = fftOps.getHE(oldId, matId/numWann, matId % numWann);
                }
            }
    return res;
}

HkFullVec WhFourierTransform::unravelFFTfullMesh(const HkFullVec &fftFullOps) const {
    assert( fftFullOps.getNumKpts() == ftp.Nf[0] * ftp.Nf[1] * ftp.Nf[2] * Nrepeat[0] * Nrepeat[1] * Nrepeat[2]);
    HkFullVec res(numWann, Nk[0] * Nk[1] * Nk[2]);
    for(unsigned a=0; a<Nk[0]; a++)
        for(unsigned b=0; b<Nk[1]; b++)
            for(unsigned c=0; c<Nk[2]; c++){
                unsigned fa = a / Nrepeat[0];
                unsigned fb = b / Nrepeat[1];
                unsigned fc = c / Nrepeat[2];
                unsigned outerId = (c % Nrepeat[2]) + Nrepeat[2] * ( (b%Nrepeat[1]) + Nrepeat[1] * (a%Nrepeat[0]));
                unsigned oldId = fc + ftp.Nf[2] * (fb + ftp.Nf[1] * (fa + ftp.Nf[0] * outerId));
                unsigned newId = c + Nk[2] * (b + Nk[1] * a);

                for(unsigned matId=0; matId<numWann*numWann; matId++){
                    res.H0(newId)[matId] = fftFullOps.getH0(oldId, matId/numWann, matId % numWann);
                    for(unsigned dir=0; dir<3; dir++){
                        res.D(newId, dir)[matId] = fftFullOps.getD(oldId, dir, matId/numWann, matId % numWann);
                        res.dH0dk(newId, dir)[matId] = fftFullOps.getdH0dk(oldId, dir, matId/numWann, matId % numWann);
                    }
                }
            }
    return res;
}
