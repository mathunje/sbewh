#ifndef SBE_WH_REQUIRED_PARAMETER_H
#define SBE_WH_REQUIRED_PARAMETER_H

#include "util/Logger.h"


enum RunMode { FermiLevel, SBE, WHtransform, Pulse };

/* - compare to Parameter
 * - the here not listed members must be parsed */
struct RequiredParameter_t {
    bool par = false;
    bool tb = false;
    bool ft = false;
    bool diag = false;
    bool pulse = false;
    bool prop = false;
    bool ksr = false;
    bool out = false;
    bool multiConfigs = false;
};

RequiredParameter_t getRequiredParameterFermiLevel();
RequiredParameter_t getRequiredParameterSBE();
RequiredParameter_t getRequiredParameterWHtransform();
RequiredParameter_t getRequiredParameterPulse();

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
