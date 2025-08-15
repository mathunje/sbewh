#include "run.h"
#include "Parameter.h"

#include "WhFourierTransform.h"
#include "util/DiagonalizatorMultiple.h"


RequiredParameter_t getRequiredParameterFermiLevel()
{
    return RequiredParameter_t {
              .par = true,
              .tb = true,
              .ft = true,
              .diag = true,
              .pulse = false,
              .prop = true,
              .ksr = false,
              .out = false,
              .multiConfigs = false
           };
}

void printEvaluation(const Parameter_t &param, const std::vector<double> &bandMin, const std::vector<double> &bandMax,
                     std::vector<double> &ev)
{
    unsigned numWann = param.tb.numWann;
    unsigned numKpts = ev.size() / numWann;
    Logger::print("Band   minimum [eV]    maximum [eV]\n");
    for(unsigned i=0; i<numWann; i++)
        Logger::print(" %3u %10lg\t%10lg\n", i, atomicUnits::to_eV(bandMin[i]), atomicUnits::to_eV(bandMax[i]));
    for(unsigned i=0; i+1<numWann; i++)
        if ( bandMax[i] < bandMin[i+1] ){
            Logger::print("Band gap of %5lg eV between band %u and %u detected\n", atomicUnits::to_eV(bandMin[i+1]-bandMax[i]), i, i+1);
            Logger::print("\tmiddle energy: %5lg eV\n", atomicUnits::to_eV(0.5 * (bandMax[i] + bandMin[i+1])) );
        }
    std::sort(ev.begin(), ev.end());
    Logger::print("Lowest eigenvalues [eV]\n");
    for(unsigned i=0; i<10; i++)
        Logger::print("  %5lg\n", atomicUnits::to_eV(ev[i]));
    Logger::print("Biggest eigenvalues [eV]\n");
    for(unsigned i=ev.size()-10; i<ev.size(); ++i)
        Logger::print("  %5lg\n", atomicUnits::to_eV(ev[i]));
    if ( param.tb.occupiedBands > 0 && param.tb.occupiedBands < numWann ){
        unsigned highestOcc = param.tb.occupiedBands * numKpts -1;
        Logger::print("For %u occupied bands:\n", param.tb.occupiedBands);
        double homo_eV = atomicUnits::to_eV(ev[highestOcc]);
        double lumo_eV = atomicUnits::to_eV(ev[highestOcc+1]);
        Logger::print("\tHOMO at %5lg eV\n", homo_eV);
        Logger::print("\tLUMO at %5lg eV\n", lumo_eV);
        Logger::print("\tFermi energy: %5lg eV\n", 0.5 * (homo_eV + lumo_eV));
    }
}


int runFermiLevelCPU(Parameter_t &param)
{
    Logger::info("Fermi level evaluation mode (CPU)\n");
    WhFourierTransform wft(param.tb, param.ft, param.prop.Nk, param.general.dim);
    GeomVector3d offset({0, 0, 0});
    const auto &shifts = wft.getFineMeshFracShifts(offset);
    param.ft.fftwPlanningFlag = FFTW_ESTIMATE;
    Logger::verbose("Planing FFT\n");
    if ( ! wft.planFFT(shifts.size()) ){
        Logger::error("Could not plan FFT\n");
        return 1;
    }
    HkVec fftMeshOps = wft.createHkVec();
    Logger::verbose("Calculating FFT\n");
    wft.transformMeshFFT(fftMeshOps, shifts, {0, 0, 0});
    Logger::verbose("Diagonalizing H0\n");

    unsigned numWann = param.tb.numWann;
    unsigned numKpts = fftMeshOps.getNumKpts();

    unsigned fmoStride = fftMeshOps.getStride();
    DiagonalizatorMultiple D(param.diag.mode, numKpts, numWann, TB_OP_CNT * fmoStride, TB_OP_POS(H0) * fmoStride );
    D.setSweepCount(param.diag.maxSweepCount);
    D.diag( fftMeshOps.data() );
    auto [bandMin, bandMax] = D.getBandLimits();
    std::vector<double> ev = D.getEigenValues();
    printEvaluation(param, bandMin, bandMax, ev);
    return 0;
}


#ifdef SBE_WH_CUDA

#include "CuWhFourierTransform.h"
#include "CuDiagonalizator.h"

void logGPUmemory()
{
    size_t free_byte;
    size_t total_byte;
    CUDA_CHECK( cudaMemGetInfo( &free_byte, &total_byte ) );
    Logger::verbose("GPU memory free: %.3lf GB | used: %.3lf GB\n", free_byte / 1024.0 / 1024.0 / 1024.0,
                     (total_byte - free_byte) / 1024.0 / 1024.0 / 1024.0);
}

int runFermiLevelGPU(Parameter_t &param)
{
    Logger::info("Fermi level evaluation mode (GPU)\n");
    cudaStream_t stream = nullptr;
    CUDA_CHECK ( cudaSetDevice(0) );
    CUDA_CHECK ( cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) );
    cuDoubleComplex * d_HkOps = nullptr;
    logGPUmemory();
    {
        Logger::verbose("Hk transform\n");
        CuWhFourierTransformHk wft(stream, param.tb, param.ft, param.prop.Noffset);
        d_HkOps = wft.createDataVec();
        wft.transform(d_HkOps, {0, 0, 0}, {0, 0, 0});
        CUDA_CHECK(cudaStreamSynchronize(stream));
        Logger::verbose("FFT calculation finished\n");
        logGPUmemory();
    }
    {
        unsigned baseCount = param.prop.Nk[0] * param.prop.Nk[1] * param.prop.Nk[2];
        CuDiagonalizator cud(stream, baseCount, param.tb.numWann, TB_OP_POS(H0) * baseCount );
        cud.diag(d_HkOps);
        std::vector<double> ev = cud.getEigenValues(true);
        auto [bandMin, bandMax] = cud.getBandLimits(ev);
        printEvaluation(param, bandMin, bandMax, ev);
        logGPUmemory();
    }
    CUDA_CHECK ( cudaFree(d_HkOps) );
    CUDA_CHECK ( cudaStreamDestroy(stream) );
    CUDA_CHECK ( cudaDeviceReset() );
    return 0;
}
#endif



int runFermiLevel(Parameter_t &param){
#ifdef SBE_WH_CUDA
    if ( param.par.useGPU ){
        return runFermiLevelGPU(param);
    } else {
        return runFermiLevelCPU(param);
    }
#else
    return runFermiLevelCPU(param);
#endif
}
