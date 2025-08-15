#ifndef SBE_WH_DIAGONALIZATION_PARAMETER_H
#define SBE_WH_DIAGONALIZATION_PARAMETER_H

#include "sbewhConfig.h"

#include "ParameterIO.h"
#include "parameter/ParameterNumeric.hpp"
#include "parameter/ParameterString.h"
#include "parallel/mpiUtil.h"
#include "util/DiagonalizatorSingle.h"
#include "util/Logger.h"

struct DiagonalizationParameter_t{
    DiagMode mode;
    unsigned maxSweepCount;
    double maxRelError;
};

class DiagonalizationParameterValidator{
    std::map<std::string, DiagMode> dmMap;
    std::string allowedModes;

    std::string diagModeStr;
    unsigned maxSweepCount;
    double maxRelError;
public:
    DiagonalizationParameterValidator(ParameterIO &paramIO);
    bool update(DiagonalizationParameter_t *p);
};

void mpi_bcast_DiagonalizationParameter(DiagonalizationParameter_t &p, const MpiParameter_t &mpi);



#endif
