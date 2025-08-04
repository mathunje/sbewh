#include "runModes.h"
#include "Parameter.h"
#include "WhFourierTransform.h"
#include "cnpy.h"


RequiredParameter_t getRequiredParameterWHtransform()
{
    return RequiredParameter_t {
              .par = true,
              .tb = true,
              .ft = true,
              .diag = false,
              .pulse = false,
              .prop = true,
              .ksr = false,
              .out = true,
              .multiConfigs = false
           };
}

std::vector<std::complex<double>> extractH0(const Parameter_t &param, unsigned Nk, const HkVec & H){
    std::vector<std::complex<double>> res;
    for(unsigned k=0; k<Nk; k++)
        for(unsigned m : param.ft.saveIndices)
            for(unsigned n : param.ft.saveIndices)
                res.push_back( H.getH0(k, m, n) );
    return res;
}

std::vector<std::complex<double>> extractD(const Parameter_t &param, unsigned Nk, const HkVec & H){
    std::vector<std::complex<double>> res;
    for(unsigned k=0; k<Nk; k++)
        for(unsigned m : param.ft.saveIndices)
            for(unsigned n : param.ft.saveIndices)
                res.push_back( H.getHE(k, m, n) - H.getH0(k, m, n) );
    return res;
}


std::vector<std::complex<double>> extractFullH0(const Parameter_t &param, unsigned Nk, const HkFullVec & H){
    std::vector<std::complex<double>> res;
    for(unsigned k=0; k<Nk; k++)
        for(unsigned m : param.ft.saveIndices)
            for(unsigned n : param.ft.saveIndices)
                res.push_back( H.getH0(k, m, n) );
    return res;
};

std::vector<std::complex<double>> extractFullD(const Parameter_t &param, unsigned Nk, const HkFullVec & H,
                                            unsigned dir){
    std::vector<std::complex<double>> res;
    for(unsigned k=0; k<Nk; k++)
        for(unsigned m : param.ft.saveIndices)
            for(unsigned n : param.ft.saveIndices)
                res.push_back( H.getD(k, dir, m, n) );
    return res;
};

std::vector<std::complex<double>> extractFulldH0dk(const Parameter_t &param, unsigned Nk, const HkFullVec & H,
                                            unsigned dir){
    std::vector<std::complex<double>> res;
    for(unsigned k=0; k<Nk; k++)
        for(unsigned m : param.ft.saveIndices)
            for(unsigned n : param.ft.saveIndices)
                res.push_back( H.getdH0dk(k, dir, m, n) );
    return res;
};

void saveWHtransform(std::string tag, const std::vector<GeomVector<double, 3>> &kMesh,
                 const std::vector< std::complex<double> > * HkOps, const Parameter_t & param, bool with_dHdk=false){
    unsigned siSize = param.ft.saveIndices.size();
    std::filesystem::path p(param.out.saveDir);
    std::string fname = std::string(p / tag) + ".npz";
    cnpy::npz_save(fname, "WanIndices", &param.ft.saveIndices[0], {siSize}, "w");
    cnpy::npz_save(fname, "lattice", unravel(param.tb.latticeVectors).data(), {3, 3}, "a");
    cnpy::npz_save(fname, "kMesh", (double*)&kMesh[0], {param.prop.Nk[0], param.prop.Nk[1], param.prop.Nk[2], 3}, "a");
    cnpy::npz_save(fname, "H0", HkOps[0].data(), {param.prop.Nk[0], param.prop.Nk[1], param.prop.Nk[2], siSize, siSize}, "a");
    std::string Ename[3] = { "Dkx", "Dky", "Dkz" };
    for(unsigned dir=0; dir<param.general.dim; dir++){
        cnpy::npz_save(fname, Ename[dir], HkOps[1+dir].data(),
                {param.prop.Nk[0], param.prop.Nk[1], param.prop.Nk[2], siSize, siSize}, "a");
    }
    if ( with_dHdk ){
        std::string Ename[3] = { "dHdkx", "dHdky", "dHdkz" };
        for(unsigned dir=0; dir<param.general.dim; dir++){
            cnpy::npz_save(fname, Ename[dir], HkOps[4+dir].data(),
                    {param.prop.Nk[0], param.prop.Nk[1], param.prop.Nk[2], siSize, siSize}, "a");
        }
    }
};


int runWHtransformCPU(Parameter_t &param){
    Logger::info("WHtransform CPU mode\n");
    WhFourierTransform wft(param.tb, param.ft, param.prop.Nk, param.general.dim);
    GeomVector3d offset({0.0, 0.0, 0.0});
    auto kMesh = wft.getKmesh(offset);
    unsigned Nk = kMesh.size();
    std::vector< std::complex<double> > HkOps[7];
    GeomVector3d E[3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} }; // test electric field

    if (param.ft.performDirectCalculation ){
        Logger::verbose("Start full Hamiltonian direct calculation\n");
        HkFullVec Hkf = wft.transformMeshDirect(offset);
        Logger::verbose("Full direct calculation finished\n");
        HkOps[0] = extractFullH0(param, Nk, Hkf);
        for(unsigned dir=0; dir<param.general.dim; dir++)
            HkOps[1+dir] = extractFullD(param, Nk, Hkf, dir);
        if ( ! param.out.disableSaving )
            saveWHtransform("trafoDirectFull", kMesh, HkOps, param);

        if ( param.ft.whTransformTestField ){
            for(unsigned dir=0; dir<param.general.dim; dir++){
                Logger::verbose("Starting direct calculation for direction %u\n", dir+1);
                HkVec Hk = wft.transformMeshDirect(offset, E[dir]);
                Logger::verbose("Direct calculation finished\n");
                if ( dir == 0 )
                    HkOps[0] = extractH0(param, Nk, Hk);
                HkOps[1+dir] = extractD(param, Nk, Hk);
            }
            if ( ! param.out.disableSaving )
                saveWHtransform("trafoDirect", kMesh, HkOps, param);
        }
    }

    const auto &shifts = wft.getFineMeshFracShifts(offset);
    Logger::verbose("Planing FFT\n");
    if ( ! wft.planFFT(shifts.size()) ){
        Logger::error("Could not plan FFT\n");
        return 1;
    }
    HkFullVec fftFullOps = wft.createHkFullVec();
    Logger::verbose("Start full FFT calculation\n");
    wft.transformMeshFFT(fftFullOps, shifts);
    Logger::verbose("Full FFT calculation finished\n");
    HkFullVec Hkf = wft.unravelFFTfullMesh(fftFullOps);
    HkOps[0] = extractFullH0(param, Nk, Hkf);
    for(unsigned dir=0; dir<param.general.dim; dir++)
        HkOps[1+dir] = extractFullD(param, Nk, Hkf, dir);
    for(unsigned dir=0; dir<param.general.dim; dir++)
        HkOps[4+dir] = extractFulldH0dk(param, Nk, Hkf, dir);
    if ( ! param.out.disableSaving )
        saveWHtransform("trafoFFTfull", kMesh, HkOps, param, true);

    if ( param.ft.whTransformTestField ) {
        HkVec fftOps = wft.createHkVec();
        for(unsigned dir=0; dir<param.general.dim; dir++){
            Logger::verbose("Starting FFT calculation for direction %u\n", dir+1);
            wft.transformMeshFFT(fftOps, shifts, E[dir]);
            Logger::verbose("FFT calculation finished\n");
            HkVec Hk = wft.unravelFFTmesh(fftOps);
            if ( dir == 0 )
                HkOps[0] = extractH0(param, Nk, Hk);
            HkOps[1+dir] = extractD(param, Nk, Hk);
        }
        if ( ! param.out.disableSaving )
               saveWHtransform("trafoFFT", kMesh, HkOps, param);
    }
    return 0;
}


#ifdef SBE_WH_CUDA

#include "CuWhFourierTransform.h"

int runWHtransformGPU(Parameter_t &param){
    Logger::info("WHtransform GPU mode\n");
    Logger::verbose("Setup cuda +  memory allocation\n");

    cudaStream_t stream = nullptr;
    CUDA_CHECK ( cudaSetDevice(0) );
    GeomVector3d offset({0.0, 0.0, 0.0});
    CUDA_CHECK ( cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) );
    {
        Logger::verbose("HkFull transform\n");
        CuWhFourierTransformHkFull wft(stream, param.tb, param.ft, param.prop.Noffset);
        cuDoubleComplex * d_HkFullOps = wft.createDataVec();
        wft.transform(d_HkFullOps, offset);
        CUDA_CHECK(cudaStreamSynchronize(stream));
        Logger::verbose("Full FFT calculation finished\n");
        HkFullVec Hkf = wft.unravelFFTmesh(d_HkFullOps);
        unsigned Nk = param.prop.Nk[0] * param.prop.Nk[1] * param.prop.Nk[2];
        std::vector< std::complex<double> > HkOps[7];
        if ( ! param.out.disableSaving ){
            HkOps[0] = extractFullH0(param, Nk, Hkf);
            for(unsigned dir=0; dir<param.general.dim; dir++){
                HkOps[1+dir] = extractFullD(param, Nk, Hkf, dir);
                HkOps[4+dir] = extractFulldH0dk(param, Nk, Hkf, dir);
            }
            auto kMesh = wft.getKmesh(offset);
            saveWHtransform("trafoCudaFull", kMesh, HkOps, param, true);
        }
    }
    if ( param.ft.whTransformTestField ){
        Logger::verbose("Hk transform\n");
        CuWhFourierTransformHk wft(stream, param.tb, param.ft, param.prop.Noffset);
        cuDoubleComplex * d_HkOps = wft.createDataVec();
        unsigned Nk = param.prop.Nk[0] * param.prop.Nk[1] * param.prop.Nk[2];
        std::array<double,3> E[3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} }; // test electric field
        std::vector< std::complex<double> > HkOps[4];
        for(unsigned dir=0; dir<param.general.dim; dir++){
            Logger::verbose("Starting FFT calculation for direction %u\n", dir+1);
            wft.transform(d_HkOps, offset, E[dir]);
            CUDA_CHECK(cudaStreamSynchronize(stream));
            Logger::verbose("FFT calculation finished\n");
            HkVec Hk = wft.unravelFFTmesh(d_HkOps);
            if ( dir == 0 )
                HkOps[0] = extractH0(param, Nk, Hk);
            HkOps[1+dir] = extractD(param, Nk, Hk);
        }
        if ( ! param.out.disableSaving ){
            auto kMesh = wft.getKmesh(offset);
            saveWHtransform("trafoCuda", kMesh, HkOps, param);
        }

        CUDA_CHECK ( cudaFree(d_HkOps) );
    }
    CUDA_CHECK ( cudaStreamDestroy(stream) );
    CUDA_CHECK ( cudaDeviceReset() );

    return 0;
}

#endif




int runWHtransform(Parameter_t &param){
#ifdef SBE_WH_CUDA
    if ( param.par.useGPU ){
        return runWHtransformGPU(param);
    } else {
        return runWHtransformCPU(param);
    }
#else
    return runWHtransformCPU(param);
#endif
}
