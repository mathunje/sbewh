#ifndef SBE_WH_PULSE_PARAMETER_H
#define SBE_WH_PULSE_PARAMETER_H

#include <stdio.h>
#include <vector>

#include "sbewhConfig.h"
#include "Logger.h"

#include "ParameterIO.h"
#include "parameter/ParameterMultiConfig.h"
#include "parameter/ParameterVector.hpp"
#include "parameter/ParameterNumeric.hpp"
#include "parameter/ParameterString.h"

#include "pulse/Pulse3D.h"
#include "pulse/Pulse1D.h"
#include "pulse/GaussianPulse1D.h"
#include "pulse/PureGauss1D.h"
#include "pulse/Sin2Pulse1D.h"
#include "pulse/Sin2RampPulse1D.h"

#include "GeomVector.hpp"

#include "unitConversion.h"
#include "mpiUtil.h"

class PulseParameterValidator;

struct PulseParameter_t {
public:
    std::string multiRunFname;
    Pulse3D pulse;
    void print(FILE *f) const;
};

class PulseParameterValidator {
private:
    const std::string &fixedInputDir;

    template<typename T> struct PulsePar_t{
        T defaultValue;
        T pulseValue[SBE_WH_MAX_PULSE_COUNT];
        const bool * setPtr[SBE_WH_MAX_PULSE_COUNT];

        bool isSet(unsigned i) const { return *setPtr[i]; }
        const T get(unsigned i) const { return ( (i>=SBE_WH_MAX_PULSE_COUNT || !(*setPtr[i]) ) ? defaultValue : pulseValue[i] ); }
    };

    PulsePar_t< std::vector<double> > pol;
    PulsePar_t< std::string> type;
    PulsePar_t<double> Emax;
    PulsePar_t<double> tCentral_fs;
    PulsePar_t<double> lambda;
    PulsePar_t<double> cep;
    PulsePar_t<double> tFWHM_fs;
    PulsePar_t<double> tStartFWHM;
    PulsePar_t<double> tEndFWHM;
    PulsePar_t<double> risingCycles;
    PulsePar_t<int> cyclesOn;

    bool addSubPulse(PulseParameter_t *p, unsigned nr, unsigned dim);

    std::map<std::string, PulseParameter_t> configurations;
    ParameterMultiConfig *configParam;
public:
    PulseParameterValidator(ParameterIO &paramIO, const std::string &fixedInputDir);

    bool newMultiConfig(std::string name);
    bool update(PulseParameter_t *p, unsigned dim=3);
    const std::map<std::string, PulseParameter_t> &getConfigurations() const { return configurations; }
};


void mpi_bcast_PulseParameter(PulseParameter_t &p, const MpiParameter_t &mp);


#endif
