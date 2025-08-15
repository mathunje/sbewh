#ifndef SBE_WH_IO_SAVE_H
#define SBE_WH_IO_SAVE_H

#include <stdio.h>
#include <filesystem>


#include "external/cnpy/cnpy.h"

#include "Parameter.h"
#include "io/util.h"
#include "git_version.h"

#include "ExpectationValues.h"

void saveInput(const Parameter_t &param);
bool saveSBEresults(const ExpectationValues_t &expValues,
                    const Parameter_t &param, const std::string &subDirectory="");

#endif
