#ifndef SBE_WH_PARAMETER_H
#define SBE_WH_PARAMETER_H

#include<string>
#include<vector>

#include "sbewhConfig.h"
#include "util/Logger.h"

#include "io/GeneralParameter.h"
#include "io/ParallelParameter.h"
#include "io/TightBindingParameter.h"
#include "io/FourierTransformParameter.h"
#include "io/DiagonalizationParameter.h"
#include "io/PulseParameter.h"
#include "io/PropagationParameter.h"
#include "io/KspaceRegionParameter.h"
#include "io/OutputParameter.h"
#include "main/runModes.h"
#include "io/util.h"
#include "parallel/MpiParameter.h"
#include "parallel/mpiUtil.h"

struct CmdInputParameter_t {
    std::vector<std::string> fnames;
    std::vector< std::string > args;
};

struct Parameter_t {
    MpiParameter_t mpi;
    CmdInputParameter_t cmdInp;
    GeneralParameter_t general;
    ParallelParameter_t par;
    TightBindingParameter_t tb;
    FourierTransformParameter_t ft;
    DiagonalizationParameter_t diag;
    PropagationParameter_t prop;
    PulseParameter_t pulse;
    KspaceRegionParameter_t ksr;
    std::map<std::string, PulseParameter_t> pulseConfigurations;
    OutputParameter_t out;
};

bool validateKdimensions(Parameter_t &param);

ParameterIO::ParseResult parseParameter(int argc, const char * argv[], Parameter_t &param, const std::string fixedInputDir="");


void mpi_bcast_Params(Parameter_t &param);


#endif
