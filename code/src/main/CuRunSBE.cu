#include "CuRunSBE.h"

double getFermiLevelGPU(cudaStream_t stream, Parameter_t &param,
                         CuWhFourierTransformHk &wft, CuDiagonalizator &cud, cuDoubleComplex * d_H, bool *ok)
{
    if ( param.tb.useCustomFermiLevel )
        return param.tb.customFermiLevel;
    Logger::verbose("Calculating Fermi level\n");
    wft.transform(d_H, {0, 0, 0}, {0, 0, 0});
    CUDA_CHECK(cudaStreamSynchronize(stream));
    unsigned nKpts = param.prop.Nk[0] * param.prop.Nk[1] * param.prop.Nk[2];
    unsigned baseOffset = nKpts * param.tb.numWann * param.tb.numWann * sizeof(std::complex<double>);
    cud.diag(d_H + TB_OP_POS(H0) * baseOffset );
    std::vector<double> ev = cud.getEigenValues(true);
    auto [bandMin, bandMax] = cud.getBandLimits(ev);
    return calcFermiLevelBase(param, bandMin, bandMax, ev, ok);
}

void propagateGPU(CuPropagator &prop, CuObserver & observer, const Parameter_t &param, GPUstate initialState)
{
    using namespace boost::numeric::odeint;
    const PropagationParameter_t & ppp = param.prop;
    double timeStep = (ppp.endTime - ppp.startTime) / (ppp.outputCount - 1);
    Logger::info("Start GPU integration\n");
    observer.reset();
    prop.reset();
    int stepCount = integrate_const(
                        make_controlled < runge_kutta_dopri5< GPUstate > >(ppp.epsAbs, ppp.epsRel),
                        std::ref(prop), initialState, ppp.startTime, ppp.endTime + timeStep/2, timeStep, std::ref(observer) );
    Logger::info("Required %d steps using GPU integration\n", stepCount);
}

int runSBE_GPU_internal(cudaStream_t stream, Parameter_t &param){
    Logger::info("SBE GPU mode\n");
    CuWhFourierTransformHk wftP(stream, param.tb, param.ft, param.prop.Noffset); // P : propagate
    CuWhFourierTransformHkFull wftF(stream, param.tb, param.ft, param.prop.Noffset); // F : full
    assert( wftP.getOpCount() <= wftF.getOpCount() );
    cuDoubleComplex * d_H = wftF.createDataVec(); // full is bigger than propagate
    unsigned baseCount = param.prop.Nk[0] * param.prop.Nk[1] * param.prop.Nk[2];
    CuDiagonalizator cud(stream, baseCount, param.tb.numWann);
    cud.setSortEig(true);
    cud.setSweepCount(param.diag.maxSweepCount);
    cud.setTolerance(param.diag.maxRelError);
    bool ok = true;
    double fermiLevel = getFermiLevelGPU(stream, param, wftP, cud, d_H, &ok);
    if ( ! ok ){
        CUDA_CHECK ( cudaFree(d_H) );
        return 1;
    }
    Logger::print("Fermi level: %.2lf eV\n", atomicUnits::to_eV(fermiLevel));
    StateContext sc(param);
    CuPropagator prop(stream, sc, wftP, d_H, cud);
    CuObserver observer(stream, sc, wftF, d_H, cud);
    GPUstate initialState = prop.calcInitialState(fermiLevel);
    if ( param.pulseConfigurations.size() ){
        for(const auto &[name, pulseParam] : param.pulseConfigurations){
            int fd = prepareMultiRun(param, name);
            if ( fd == -1)
                continue;
            Logger::info("Configuration: %s\n", name.c_str());
            param.pulse = pulseParam;
            propagateGPU(prop, observer, param, initialState);
            if ( ! param.out.disableSaving ){
                const ExpectationValues_t & expValues = observer.collectExpectationValues();
                if ( ! saveSBEresults(expValues, param, name) ){
                    CUDA_CHECK ( cudaFree(d_H) );
                    return 1;
                }
            }
            finishMultiRun(param, fd);
        }
    } else {
        propagateGPU(prop, observer, param, initialState);
        if ( ! param.out.disableSaving ){
            const ExpectationValues_t & expValues = observer.collectExpectationValues();
            if ( ! saveSBEresults(expValues, param) ){
                CUDA_CHECK ( cudaFree(d_H) );
                return 1;
            }
        }
    }
    CUDA_CHECK ( cudaFree(d_H) );
    return 0;
}

int runSBE_GPU(Parameter_t &param)
{
    cudaStream_t stream = nullptr;
    CUDA_CHECK ( cudaSetDevice(0) );
    CUDA_CHECK ( cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) );
    int res = runSBE_GPU_internal(stream, param);
    CUDA_CHECK ( cudaStreamDestroy(stream) );
    CUDA_CHECK ( cudaDeviceReset() );
    return res;
}
