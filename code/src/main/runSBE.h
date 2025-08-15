#ifndef SBE_WH_RUN_SBE_H
#define SBE_WH_RUN_SBE_H

#include "sbewhConfig.h"
#include "Parameter.h"
#include "WhFourierTransform.h"
#include "util/DiagonalizatorMultiple.h"
#include "PropagatorCPU.h"
#include "ObserverCPU.h"

#include "io/save.h"
#include "parallel/mpiUtil.h"

#include<boost/numeric/odeint.hpp>
#include <cassert>
#include<thread>

#include <sys/file.h>
#include <unistd.h>
#include <errno.h>

#ifdef SBE_WH_CUDA
#include "CuRunSBE.h"
#endif

RequiredParameter_t getRequiredParameterSBE();

int prepareMultiRun(const Parameter_t &param, const std::string &name);
void finishMultiRun(const Parameter_t &param, int fd);
bool assignProcesses(Parameter_t &param, unsigned &NrStart, unsigned &NrCount);
double calcFermiLevelBase(Parameter_t &param, std::vector<double> &bandMin, std::vector<double> &bandMax,
                          std::vector<double> ev, bool *ok);
int runSBE(Parameter_t &param);

#endif
