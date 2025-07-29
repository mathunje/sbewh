#ifndef SBE_WH_OUTPUT_PARAMETER_H
#define SBE_WH_OUTPUT_PARAMETER_H

#include<stdio.h>

#include<string>
#include<filesystem>

#include "sbewhConfig.h"
#include "Logger.h"

#include "ParameterIO.h"
#include "parameter/ParameterBool.h"
#include "parameter/ParameterString.h"

#include "mpiUtil.h"

struct OutputParameter_t {
    bool disableSaving;
    bool saveAsNpz;
    bool saveInput;
    bool saveWannierInput;
    std::string baseDir;
    std::string run;
    std::string saveDir; // baseDir/run if run set otherwise baseDir/runNumber
};


class OutputParameterValidator {
private:
    bool disableSaving;
    bool saveAsNpz;
    bool saveInput;
    bool saveWannierInput;
    const bool * saveWannierInputSetPtr;
    std::string baseDir;
    std::string run;
public:
    OutputParameterValidator(ParameterIO &paramIO);
    bool update(OutputParameter_t *p, bool selectLastRun);
};


void mpi_bcast_OutputParameter(OutputParameter_t &p, const MpiParameter_t &mp);


#endif
