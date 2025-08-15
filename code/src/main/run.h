#ifndef SBE_WH_RUN_H
#define SBE_WH_RUN_H

#include "sbewhConfig.h"
#include "runModes.h"
#include "util/Logger.h"
#include "runSBE.h"

int runFermiLevel(Parameter_t &param);
int runPulse(Parameter_t &param);
int runWHtransform(Parameter_t &param);

#endif
