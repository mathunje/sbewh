#include "runSBE.h"

RequiredParameter_t getRequiredParameterSBE()
{
    return RequiredParameter_t {
              .par = true,
              .tb = true,
              .ft = true,
              .diag = true,
              .pulse = true,
              .prop = true,
              .ksr = true,
              .out = true,
              .multiConfigs = true
           };
}


int prepareMultiRun(const Parameter_t &param, const std::string &name)
{
    int retVal = 0;
    if ( param.mpi.rank == param.mpi.root ){
        std::string lockFile = std::filesystem::path(param.out.saveDir) / ("." + name);
        int fd = open(lockFile.c_str(), O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
        if ( errno == EEXIST ){
            int fd = open(lockFile.c_str(), O_WRONLY);
            if ( fd == -1 ){
                retVal = -1;
                Logger::error("Could not open lock file '%s'\n", lockFile.c_str());
            } else {
                retVal = fd;
                if ( flock(fd, LOCK_EX | LOCK_NB) && errno == EWOULDBLOCK ){
                    close(fd);
                    retVal = -1;
                } else {
                    off_t fsize = lseek(fd, 0, SEEK_END);
                    if ( fsize ){
                        close(fd);
                        retVal = -1;
                    }
                }
            }
        } else {
            if ( fd == -1 ){
                Logger::error("Could not create lock file '%s'\n", lockFile.c_str() );
                retVal = -1;
            } else {
                retVal = fd;
                if ( flock(fd, LOCK_EX | LOCK_NB) && errno == EWOULDBLOCK ){
                    printf("Locked\n");
                    close(fd);
                    retVal = -1;
                }
            }
        }
    }
    mpi_bcast_trivial(retVal, param.mpi);
    return retVal;
}

void finishMultiRun(const Parameter_t &param, int fd)
{
    if ( param.mpi.rank == param.mpi.root){
        if ( fd != -1 ){
            char s = '1';
            if ( write(fd, &s, 1) != 1 )
                Logger::warn("Finished flag could not be written\n");
            flock(fd, LOCK_UN);
            close(fd);
        }
    }
}

bool assignProcesses(Parameter_t &param, unsigned &NrStart, unsigned &NrCount){
    const std::array<unsigned, 3> &Nk = param.prop.Nk;
    const std::array<unsigned, 3> &Nf = param.ft.Nf;
    unsigned Nrepeat = (Nk[0] / Nf[0]) * (Nk[1] / Nf[1]) * (Nk[2] / Nf[2]);
    unsigned Nprocess = ( Nrepeat + param.mpi.size - 1) / param.mpi.size;
    NrStart = param.mpi.rank * Nprocess;
    unsigned NrEnd = std::min( (param.mpi.rank+1) * Nprocess, Nrepeat);
    // avoid having processes without work:
    if ( (param.mpi.size -1) * Nprocess >= Nrepeat){
        if ( Nprocess == 1 ){
            if ( param.mpi.rank == param.mpi.root )
                Logger::error( "More processes than number of shifted grids\n");
            return false;
        }
        Nprocess -= 1;
        NrStart = param.mpi.rank * Nprocess;
        NrEnd = (param.mpi.rank + 1 == param.mpi.size) ? Nrepeat : (param.mpi.rank+1) * Nprocess;
    }
    NrCount = NrEnd - NrStart;
    return true;
}


double calcFermiLevelBase(Parameter_t &param, std::vector<double> &bandMin, std::vector<double> &bandMax,
                          std::vector<double> ev, bool *ok)
{
    *ok = true;
    unsigned numWann = param.tb.numWann;
    unsigned numKpts = ev.size() / numWann;
    std::sort(ev.begin(), ev.end());
    if ( param.tb.occupiedBelowBandGap ){
        double fermiLevel = 0;
        unsigned gapCount = 0;
        for(unsigned i=0; i+1<numWann; i++)
            if ( bandMax[i] + 3 * param.tb.gsTemp  < bandMin[i+1] ){
                gapCount++;
                fermiLevel = 0.5 * (bandMax[i] + bandMin[i+1]);
            }
        if ( gapCount != 1 ){
            Logger::error("Detected %u != 1 band gaps\n", gapCount);
            *ok = false;
        }
        return fermiLevel;
    }
    assert( param.tb.occupiedBands > 0 && param.tb.occupiedBands < numWann );
    unsigned highestOcc = param.tb.occupiedBands * numKpts -1;
    double guessHomo = ev[highestOcc];
    double guessLumo = ev[highestOcc+1];
    double fermiLevelGuess = 0.5 * (guessHomo + guessLumo);
    if ( guessHomo + 3 * param.tb.gsTemp > guessLumo){
        Logger::error("Fermi level calculation for metals is currently not implemented\n");
        *ok = false;
    }
    return fermiLevelGuess;
}

double getFermiLevelCPU(Parameter_t &param, bool *ok)
{
    if ( param.tb.useCustomFermiLevel )
        return param.tb.customFermiLevel;
    Logger::verbose("Calculating Fermi level\n");
    WhFourierTransform wft(param.tb, param.ft, param.prop.Nk, param.general.dim);
    GeomVector3d offset({0, 0, 0});
    const auto &shifts = wft.getFineMeshFracShifts(offset);
    if ( ! wft.planFFT(shifts.size(), FFTW_ESTIMATE) ){
        Logger::error("calcFermiLevel: Could not initalize FFT\n");
        *ok = false;
        return 0;
    }
    HkVec fftMeshOps = wft.createHkVec();
    std::array<double,3> Ezero({0, 0, 0});
    wft.transformMeshFFT(fftMeshOps, shifts, Ezero);
    unsigned numWann = param.tb.numWann;
    unsigned fmoStride = fftMeshOps.getStride();
    unsigned numKpts = fftMeshOps.getNumKpts();
    DiagonalizatorMultiple D(param.diag.mode, numKpts, numWann, TB_OP_CNT * fmoStride, TB_OP_POS(H0) * fmoStride );
    D.setSweepCount(param.diag.maxSweepCount);
    D.diag(fftMeshOps.data());
    auto [bandMin, bandMax] = D.getBandLimits();
    std::vector<double> ev = D.getEigenValues();
    return calcFermiLevelBase(param, bandMin, bandMax, ev, ok);
}

bool saveCPUresultsIfDesired(const Parameter_t &param, const ObserverCPU &observer, std::string configName="")
{
    if ( param.out.disableSaving )
        return true;
    bool ok = true;
    if ( param.mpi.rank == param.mpi.root ){
        const ExpectationValues_t & expValues = observer.collectExpectationValues();
        ok = saveSBEresults(expValues, param, configName);
    }
    mpi_bcast_trivial(ok, param.mpi);
    return ok;
}

void propagateCPUintertwined(PropagatorCPU &prop, ObserverCPU & observer, const Parameter_t &param, CPUstate initialState)
{
    using namespace boost::numeric::odeint;
    const PropagationParameter_t & ppp = param.prop;
    double timeStep = (ppp.endTime - ppp.startTime) / (ppp.outputCount - 1);
    Logger::info("Start CPU integration\n");
    prop.reset();
    observer.reset();
    int stepCount = integrate_const(
                        make_controlled < runge_kutta_dopri5< CPUstate > >(ppp.epsAbs, ppp.epsRel),
                        std::ref(prop), initialState, ppp.startTime, ppp.endTime + timeStep/2, timeStep, std::ref(observer) );
    Logger::info("Required %d steps using CPU integration\n", stepCount);
    observer.normalizeExpectationValues();
}

int intertwinedIntegrations(Parameter_t &param, double fermiLevel)
{
    unsigned NrProcStart, NrProcCount;
    if ( ! assignProcesses(param, NrProcStart, NrProcCount) )
        return 1;
    StateContext sc(param, NrProcStart, NrProcCount);
    WhFourierTransform wft(param.tb, param.ft, param.prop.Nk, param.general.dim);
    Logger::verbose("Planing FFT for synchronous integration\n");
    int failCount = 0;
    if ( ! wft.planFFT(sc.NrProcCount) )
        failCount = 1;
    mpi_sum(failCount, param.mpi);
    if ( failCount ){
        if ( param.mpi.rank == param.mpi.root)
            Logger::error("Could not plan FFT for %d processes\n", failCount);
        return 1;
    }
    PropagatorCPU prop(sc, wft);
    ObserverCPU observer(sc, wft);
    std::vector<std::complex<double>> initialState = prop.calcInitialState(fermiLevel);
    if ( param.pulseConfigurations.size() ){
        for(const auto &[name, pulseParam] : param.pulseConfigurations){
            int fd = prepareMultiRun(param, name);
            if ( fd == -1 )
                continue;
            Logger::info("Configuration: %s\n", name.c_str());
            param.pulse = pulseParam;
            propagateCPUintertwined(prop, observer, param, initialState);
            if ( ! saveCPUresultsIfDesired(param, observer, name) )
                return 1;
            finishMultiRun(param, fd);
        }
    } else {
        propagateCPUintertwined(prop, observer, param, initialState);
        if ( ! saveCPUresultsIfDesired(param, observer) )
            return 1;
    }
    return 0;
}

#ifdef SBE_WH_MPI
void assignTasks(double sleepTime, int mpiRoot, int mpiSize, int Nr, bool *rootFinished, bool *nValid, int *nNew)
{
    std::vector<int> intRound(mpiSize, 0);
    std::vector<MPI_Request> request(mpiSize);
    unsigned dummy;
    for(int i=0; i<mpiSize; i++)
        if ( i != mpiRoot )
            MPI_Irecv(&dummy, 1, MPI_INT, i, intRound[i], MPI_COMM_WORLD, &request[i]);
    bool allDone = false;
    int nextN = mpiSize;
    while ( !allDone ){
        sleep(sleepTime);
        allDone = true;
        for(int i=0; i<mpiSize; i++){
            if (intRound[i] < 0)
                continue;
            allDone = false;
            if ( i == mpiRoot ){
                if ( *rootFinished ){
                    *rootFinished = false;
                    *nNew = nextN;
                    *nValid = true;
                    if ( nextN >= Nr )
                        intRound[i] = -1;
                    nextN++;
                }
            } else {
                int finished;
                MPI_Test(&request[i], &finished, MPI_STATUS_IGNORE);
                if ( finished ){
                    MPI_Send(&nextN, 1, MPI_INT, i, intRound[i], MPI_COMM_WORLD);
                    if ( nextN < Nr ){
                        intRound[i]++;
                        MPI_Irecv(&dummy, 1, MPI_INT, i, intRound[i], MPI_COMM_WORLD, &request[i]);
                    } else {
                        intRound[i] = -1;
                    }
                    nextN++;
                }
            }
        }
    }
}

void propagateCPUindependent(StateContext &sc, PropagatorCPU &prop, ObserverCPU & observer, const Parameter_t &param, double fermiLevel)
{
    const std::array<unsigned, 3> &Nk = param.prop.Nk;
    const std::array<unsigned, 3> &Nf = param.ft.Nf;
    unsigned Nrepeat = int( (Nk[0] / Nf[0]) * (Nk[1] / Nf[1]) * (Nk[2] / Nf[2]) );
    using namespace boost::numeric::odeint;
    const PropagationParameter_t & ppp = param.prop;
    double timeStep = (ppp.endTime - ppp.startTime) / (ppp.outputCount - 1);
    Logger::info("Start CPU integration\n");
    prop.reset();
    observer.reset();
    int Nr = param.mpi.rank;
    int intRound = 0;
    bool calcFinished = false;
    bool nValid = false;
    int NrNext = 0;
    std::thread thread_assign;
    if ( param.mpi.rank == param.mpi.root )
        thread_assign = std::thread(assignTasks, param.par.rmaPollTime, param.mpi.root, param.mpi.size, Nrepeat,
                                    &calcFinished, &nValid, &NrNext);
    MPI_Request request;
    while ( Nr < (int)Nrepeat ){
        // actual calculation
        sc.NrProcStart = Nr; // stateContext is used by Propagator and Observer
        std::vector<std::complex<double>> initialState = prop.calcInitialState(fermiLevel);
        int stepCount = integrate_const(
                        make_controlled < runge_kutta_dopri5< CPUstate > >(ppp.epsAbs, ppp.epsRel),
                        std::ref(prop), initialState, ppp.startTime, ppp.endTime + timeStep/2, timeStep, std::ref(observer) );
        observer.setLocalMerging(true);
        if ( param.general.verbose ){
            unsigned Nmod = std::max(Nrepeat/20, 1u);
            if ( param.general.veryVerbose || Nr % Nmod == 0 ){
                Logger::print("Finished %d/%d\n", Nr, Nrepeat);
            }
        }
        //get net Nr
        if ( param.mpi.rank == param.mpi.root ){
            calcFinished = true;
            while( ! nValid )
                sleep(param.par.rmaPollTime);
            Nr = NrNext;
            nValid = false;
        } else {
            int one = 1;
            MPI_Isend(&one, 1, MPI_INT, param.mpi.root, intRound, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
            MPI_Recv(&Nr, 1, MPI_INT, param.mpi.root, intRound, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            intRound++;
        }
    }
    if ( param.mpi.rank == param.mpi.root )
        thread_assign.join();

    observer.mpi_reduceExpectationValues();
    observer.normalizeExpectationValues();
}


int independentIntegrations(Parameter_t &param, double fermiLevel)
{
    const std::array<unsigned, 3> &Nk = param.prop.Nk;
    const std::array<unsigned, 3> &Nf = param.ft.Nf;
    unsigned Nrepeat = (Nk[0] / Nf[0]) * (Nk[1] / Nf[1]) * (Nk[2] / Nf[2]);
    if ( param.mpi.size > (int)Nrepeat ){
        Logger::error("More processes (%d) than independent k-shifts (%u)\n", param.mpi.size, Nrepeat);
        return 1;
    }
    StateContext sc(param, 0, 1);
    WhFourierTransform wft(param.tb, param.ft, param.prop.Nk, param.general.dim);
    Logger::verbose("Planing FFT\n");
    if ( ! wft.planFFT( 1 ) ){
        Logger::error("Could not plan FFT\n");
        return 1;
    }
    PropagatorCPU prop(sc, wft);
    ObserverCPU observer(sc, wft);
    if ( param.pulseConfigurations.size() ){
        for(const auto &[name, pulseParam] : param.pulseConfigurations){
            int fd = prepareMultiRun(param, name);
            if ( fd == -1)
                continue;
            Logger::info("Configuration: %s\n", name.c_str());
            param.pulse = pulseParam;
            propagateCPUindependent(sc, prop, observer, param, fermiLevel);
            if ( ! saveCPUresultsIfDesired(param, observer, name) )
                return 1;
            finishMultiRun(param, fd);
        }
    } else {
        propagateCPUindependent(sc, prop, observer, param, fermiLevel);
        if ( ! saveCPUresultsIfDesired(param, observer) )
            return 1;
    }
    return 0;
}
#endif

int runSBE_CPU(Parameter_t &param){
    double fermiLevel = 0;
    bool ok = true;
    if ( param.mpi.rank == param.mpi.root ){
        Logger::info("SBE CPU mode\n");
        fermiLevel = getFermiLevelCPU(param, &ok);
        if ( ok )
            Logger::print("Fermi level: %.2lf eV\n", atomicUnits::to_eV(fermiLevel));
    }
    mpi_bcast_trivial(ok, param.mpi);
    mpi_bcast_trivial(fermiLevel, param.mpi);
    if ( ! ok )
        return 1;
#ifdef SBE_WH_MPI
    if ( param.par.independentMpiIntegrations ){
        return independentIntegrations(param, fermiLevel);
    } else {
        return intertwinedIntegrations(param, fermiLevel);
    }
#else
    return intertwinedIntegrations(param, fermiLevel);
#endif
}

int runSBE(Parameter_t &param){
#ifdef SBE_WH_CUDA
    if ( param.par.useGPU )
        return runSBE_GPU(param);
#endif
    return runSBE_CPU(param);
}
