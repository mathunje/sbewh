#ifndef SBE_WH_TIGHTBINDING_PARAMETER_H
#define SBE_WH_TIGHTBINDING_PARAMETER_H

#include <stdio.h>
#include <math.h>
#include <filesystem>

#include <vector>
#include <map>
#include <set>
#include <complex>

#include "sbewhConfig.h"
#include "ParameterIO.h"
#include "parameter/ParameterNumeric.hpp"
#include "parameter/ParameterVector.hpp"
#include "parameter/ParameterString.h"
#include "parameter/ParameterBool.h"
#include "util/Logger.h"
#include "util/w90util.h"

#include "util/GeomVector.hpp"
#include "util/unitConversion.h"

#include "parallel/mpiUtil.h"

struct TightBindingParameter_t {
    unsigned numWann;
    double gsTemp;
    unsigned occupiedBands;
    bool occupiedBelowBandGap;
    double customFermiLevel;
    bool useCustomFermiLevel;

    std::array<GeomVector3d, 3> latticeVectors;
    std::unordered_map<CellIndex, W90_realSpaceOperators> realSpaceOps;
    std::string wannierTbFname;
    std::string wannierWsvecFname;
    std::array<unsigned, 3> Nw;
    CellIndex minCellIndices;

    void print(FILE *f, unsigned dim) const;
};

class TightBindingParameterValidator{
private:
    const std::string &fixedInputDir;
    std::array<GeomVector3d, 3> latticeVectors;
    std::string wannierSeed;
    bool symmetrizeHamiltonian;
    bool onlyRealMatrixElements;
    unsigned occupiedBands;
    bool occupiedBelowBandGap;
    double gsTemp;
    double customFermiLevel;
    const bool *customFermiLevelSetPtr;
    bool evalWannierFiles(TightBindingParameter_t * p, unsigned dim);
public:
    TightBindingParameterValidator(ParameterIO &paramIO, const std::string &fixedInputDir);
    bool update(TightBindingParameter_t *p, unsigned dim);
};

void mpi_bcast_TightBindingParameter(TightBindingParameter_t &p, const MpiParameter_t &mp);

#endif
