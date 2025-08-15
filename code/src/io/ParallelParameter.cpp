#include "io/ParallelParameter.h"

ParallelParameterValidator::ParallelParameterValidator(ParameterIO &paramIO)
{
    paramIO.startGroup("Parallelization", "for CPU, Parallelization using indepdendent MPI integrations and no openmp threads is probably the fastest option\n");
    paramIO.registerParam(new ParameterNumeric<unsigned>("maxThreadCount", "maximum number of threads for omp", maxThreadCount, 1000000));
    paramIO.registerParam(new ParameterNumeric<double>("rmaPollTime", "wait time [s] until next mpi sync during independent MPI integrations is performed", rmaPollTime, 0.1));
    paramIO.registerParam(new ParameterBool("independentMpiIntegrations", "if enabled the integrations for different k points are "
                                            "performed seperately and the results are finally merged", independentMpiIntegrations, true));
    paramIO.registerParam(new ParameterBool("useGPU", "if enabled all calculations are performed with cuda", useGPU, false));
    paramIO.finishGroup();
}

bool ParallelParameterValidator::update(ParallelParameter_t *p, const MpiParameter_t &mpi)
{
    p->maxThreadCount = maxThreadCount;
#ifndef SBE_WH_CUDA
    if ( useGPU ){
        Logger::warn("useGPU is set to false as application was compiled without cuda support\n");
        useGPU = false;
    }
#endif
    p->rmaPollTime = rmaPollTime;
    p->independentMpiIntegrations = independentMpiIntegrations;
    p->useGPU = useGPU;
    if ( p->independentMpiIntegrations && p->useGPU ){
        Logger::warn("Disabling independent MPI integrations as they are currently not implemented for GPU\n");
        p->independentMpiIntegrations = false;
    }
#ifdef SBE_WH_MPI
    if ( p->independentMpiIntegrations ){
        if ( mpi.threadModeProvided < MPI_THREAD_SERIALIZED ){
            Logger::error("Given MPI implementation does not provide MP_SERIALIZED, which is required for independent mpi integrations\n");
            return false;
        }
    }
#endif
    return true;
}

void mpi_bcast_ParallelParameter(ParallelParameter_t &p, const MpiParameter_t &mp)
{
    mpi_bcast_trivial(p.maxThreadCount, mp);
    mpi_bcast_trivial(p.rmaPollTime, mp);
    mpi_bcast_trivial(p.independentMpiIntegrations, mp);
    mpi_bcast_trivial(p.useGPU, mp);
}
