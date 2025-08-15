#ifndef SBE_WH_PARALLEL_PARAMETER_H
#define SBE_WH_PARALLEL_PARAMETER_H

#include "sbewhConfig.h"
#include "util/Logger.h"
#include "ParameterIO.h"
#include "parameter/ParameterNumeric.hpp"
#include "parameter/ParameterBool.h"

#include "parallel/MpiParameter.h"
#include "parallel/mpiUtil.h"

struct ParallelParameter_t{
    unsigned maxThreadCount;
    double rmaPollTime;
    bool independentMpiIntegrations;
    bool useGPU;
};

class ParallelParameterValidator{
private:
    unsigned maxThreadCount;
    double rmaPollTime;
    bool independentMpiIntegrations;
    bool useGPU;
public:
    ParallelParameterValidator(ParameterIO &paramIO);
    bool update(ParallelParameter_t *p, const MpiParameter_t &mpi);
};

void mpi_bcast_ParallelParameter(ParallelParameter_t &p, const MpiParameter_t &mp);

#endif
