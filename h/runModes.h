#ifndef SBE_WH_RUN_MODES_H
#define SBE_WH_RUN_MODES_H

#include "sbewhConfig.h"
#include "Logger.h"
#include "RequiredParameter.h"


enum RunMode { FermiLevel, SBE, WHtransform, Pulse };

struct Parameter_t;

RequiredParameter_t getRequiredParameterFermiLevel();
int runFermiLevel(Parameter_t &param);

RequiredParameter_t getRequiredParameterSBE();
int runSBE(Parameter_t &param);

int runWHtransform(Parameter_t &param);
RequiredParameter_t getRequiredParameterWHtransform();


RequiredParameter_t getRequiredParameterPulse();
int runPulse(Parameter_t &param);

inline RequiredParameter_t getRequiredParameter(RunMode runMode){
    switch(runMode){
        case FermiLevel: return getRequiredParameterFermiLevel();
        case SBE: return getRequiredParameterSBE();
        case WHtransform: return getRequiredParameterWHtransform();
        case Pulse: return getRequiredParameterPulse();
        default:
            Logger::warn("This warning should never occur: Please update runModes.h");
            return RequiredParameter_t();
    }
}

#endif
