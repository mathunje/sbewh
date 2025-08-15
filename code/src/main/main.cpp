#include "sbewhConfig.h"
#include "git_version.h"
#include "util/Logger.h"
#include "run.h"

#include "Parameter.h"
#include "io/save.h"
#include "io/util.h"

#include "parallel/mpiUtil.h"


#ifdef SBE_WH_OPENMP
#include <omp.h>
#include <thread>
#endif


int main(int argc, const char *argv[])
{
    Parameter_t param;
    char ** ca = const_cast<char**>(argv);
    mpi_init(param.mpi, argc, ca);
    ParameterIO::ParseResult parseResult = ParameterIO::ParseResult::Ok;
    if (param.mpi.rank == param.mpi.root){
        Logger::setVerbosity(Logger::Severity::INFO);
        parseResult = parseParameter(argc, argv, param);
    }
    mpi_bcast_trivial(parseResult, param.mpi);
    if ( parseResult==ParameterIO::ParseResult::Fail ||
         parseResult==ParameterIO::ParseResult::HelpRequested ){
        mpi_finalize();
        return (parseResult != ParameterIO::ParseResult::HelpRequested);
    }
    if (param.mpi.rank == param.mpi.root){
        char name[50];
        sprintf(name, "SBE in Wannier Houston basis (version %u.%u)", SBE_WH_VERSION_MAJOR, SBE_WH_VERSION_MINOR);
        Logger::helloWorld(name, kGitHash);
    }
    mpi_bcast_Params(param);
    if ( param.general.dryRun ){
        Logger::info("Dry run (parameter parsing) finished\n");
        mpi_finalize();
        return 0;
    }
    if ( param.mpi.size > 1)
        Logger::verbose("MPI size: %d\n", param.mpi.size);
    RequiredParameter_t parsedParams = getRequiredParameter(param.general.runMode);
    if ( parsedParams.par ){
#ifdef SBE_WH_OPENMP
        unsigned nthreads = std::thread::hardware_concurrency();
        omp_set_nested(0);
        Logger::verbose("Number of hardware threads (on root): %u\n", nthreads);
        if ( param.par.maxThreadCount < nthreads ){
            Logger::info("Limiting calculation to %u threads\n", param.par.maxThreadCount);
            omp_set_num_threads(param.par.maxThreadCount);
        }
        fftw_init_threads();
        fftw_plan_with_nthreads(omp_get_num_threads());
#endif
    }
    if ( parsedParams.out ){
        if ( ! param.out.disableSaving && ! param.general.finishMultiRuns ){
            if ( param.mpi.rank == param.mpi.root ){
                createOutputDirectory(param.out);
                if ( param.out.saveInput )
                    saveInput(param);
            }
        }
    }
    int res = 0;
    if ( param.mpi.size > 1 ){
        switch (param.general.runMode){
            case WHtransform:
            case FermiLevel:
                if ( param.mpi.rank == param.mpi.root )
                    Logger::error("MPI mode is only for run mode SBE on CPU implemented\n");
                break;
            case SBE:
                if ( param.par.useGPU ){
                    if ( param.mpi.rank == param.mpi.root )
                        Logger::error("MPI mode is only for run mode SBE on CPU implemented\n");
                } else {
                    res = runSBE(param);
                }
                break;
            case Pulse:
                if ( param.mpi.rank == param.mpi.root )
                    res = runPulse(param);
                break;
            default:
                if ( param.mpi.rank == param.mpi.root )
                    Logger::error("Unexpected run mode\n");
        }
    } else {
        switch (param.general.runMode){
            case WHtransform:
                res = runWHtransform(param);
                break;
            case FermiLevel:
                res = runFermiLevel(param);
                break;
            case SBE:
                res = runSBE(param);
                break;
            case Pulse:
                res = runPulse(param);
                break;
            default:
                Logger::error("Unexpected run mode\n");
        }
    }

#ifdef SBE_WH_OPENMP
    fftw_cleanup_threads();
#endif
    mpi_finalize();
    if ( param.mpi.rank == param.mpi.root )
        Logger::info("Program execution finished\n");
    return res;
}
