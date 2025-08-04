#ifndef SBE_WH_FFT_PARAMETER_H
#define SBE_WH_FFT_PARAMETER_H

#include "sbewhConfig.h"
#include "Logger.h"
#include "ParameterIO.h"
#include "parameter/ParameterVector.hpp"
#include "parameter/ParameterNumeric.hpp"
#include "parameter/ParameterBool.h"
#include "parameter/ParameterString.h"

#include "io/TightBindingParameter.h"
#include "mpiUtil.h"

#include <fftw3.h>
#include <map>

struct FourierTransformParameter_t{
    bool performDirectCalculation;
    bool whTransformTestField;;
    std::vector<unsigned> saveIndices;
    std::array<unsigned, 3> Nf;
    unsigned maxNfPrime;
    unsigned fftwPlanningFlag;
    double planningTimeLimit;
    void print(FILE *f, unsigned dim) const;
};

class FourierTransformParameterValidator{
private:
    bool performDirectCalculation;
    bool whTransformTestField;
    std::vector<unsigned> saveIndices;
    std::vector<unsigned> Nf;
    unsigned maxNfPrime;
    std::string fftPlanningMode;
    double planningTimeLimit;
    static std::map<std::string, unsigned> planningModeMap;
public:
    FourierTransformParameterValidator(ParameterIO &paramIO);
    bool update(FourierTransformParameter_t *p, const TightBindingParameter_t &tb, unsigned dim);

    std::string static planningFlagToStr(unsigned flag);
};

void mpi_bcast_FourierTransformParameter(FourierTransformParameter_t &p, const MpiParameter_t &mp);

#endif
