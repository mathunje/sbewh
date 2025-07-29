#ifndef SBE_WH_GENERAL_PARAMETER_H
#define SBE_WH_GENERAL_PARAMETER_H

#include "sbewhConfig.h"
#include "Logger.h"
#include "ParameterIO.h"
#include "parameter/ParameterNumeric.hpp"
#include "parameter/ParameterString.h"
#include "parameter/ParameterBool.h"
#include "mpiUtil.h"

#include "runModes.h"

struct GeneralParameter_t{
    bool dryRun;
    unsigned dim;
    bool verbose;
    bool veryVerbose;
    bool timeIt;
    RunMode runMode;
    bool finishMultiRuns;
};

class GeneralParameterValidator{
private:
    std::map<std::string, RunMode> rmMap;
    bool dryRun;
    unsigned dim;
    bool verbose;
    bool veryVerbose;
    bool timeIt;
    std::string runModeStr;
    bool finishMultiRuns;
public:
    GeneralParameterValidator(ParameterIO &paramIO);
    bool update(GeneralParameter_t *p);
};

void mpi_bcast_GeneralParameter(GeneralParameter_t &p, const MpiParameter_t &mp);

#endif
