#ifndef SBE_WH_KSPACE_REGION_PARAMETR_H
#define SBE_WH_KSPACE_REGION_PARAMETR_H

#include "sbewhConfig.h"
#include "Logger.h"
#include "ParameterIO.h"
#include "parameter/ParameterBool.h"
#include "parameter/ParameterNumeric.hpp"

#include "GeomVector.hpp"

#include "mpiUtil.h"

#include<vector>
#include<array>
#include<string>


struct KspaceRegion_t {
    std::string name;
    bool movingFrame;
    double fwhm;
    std::array<unsigned, 3> sampleDim;
    GeomVector3d kFrac;
    GeomVector3d dkFrac[3];
};

struct KspaceRegionParameter_t {
    bool storeDensityMatrix;
    unsigned densityOutputStride;
    bool storeExpectationValues;
    std::vector<KspaceRegion_t> regions;
    void print(FILE *f) const;
};

class KspaceRegionParameterValidator{
private:
    std::map<std::string, std::string> params;
    bool storeDensityMatrix;
    unsigned densityOutputStride;
public:
    KspaceRegionParameterValidator(ParameterIO &paramIO);
    bool update(KspaceRegionParameter_t *p);
};


void mpi_bcast_KspaceRegionParameter(KspaceRegionParameter_t &p, const MpiParameter_t &mp);

#endif
