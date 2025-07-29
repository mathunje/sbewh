#ifndef SBE_WH_PROPAGATION_PARAMETER_H
#define SBE_WH_PROPAGATION_PARAMETER_H

#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <numbers>

#include "sbewhConfig.h"
#include "Logger.h"
#include "ParameterIO.h"
#include "parameter/ParameterVector.hpp"
#include "parameter/ParameterNumeric.hpp"
#include "parameter/ParameterBool.h"
#include "parameter/ParameterString.h"

#include "unitConversion.h"

#include "io/PulseParameter.h"
#include "mpiUtil.h"

struct PropagationParameter_t {
    unsigned outputCount;
    double startTime;
    double endTime;
    double epsAbs;
    double epsRel;

    double relaxationF2;
    double soothingWidth;
    double occupationSmearingWidth;

    std::array<unsigned, 3> Nk;
    std::array<unsigned, 3> Noffset; // Noffset[d] * ft.Nf[d] = N[k]
    GeomVector3d kGlobalShift;
    bool allowNkAdaption;

    void print(FILE *f);
};


class PropagationParameterValidator{
    double T2_fs;
    double soothingWidth_eV;
    double occupationSmearingWidth_eV;
    const bool *occupationSmearingWidthSetPtr;
    std::vector<unsigned> Nk;
    std::vector<double> kGlobalShift;
    bool allowNkAdaption;
    unsigned outputCount;
    double epsAbs;
    double epsRel;
public:
    PropagationParameterValidator(ParameterIO &paramIO);
    bool update(PropagationParameter_t *p, const PulseParameter_t &pulse, unsigned dim);
};

void mpi_bcast_PropagationParameter(PropagationParameter_t &p, const MpiParameter_t &mp);

#endif
