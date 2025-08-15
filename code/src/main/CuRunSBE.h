#ifndef SBE_WH_CU_RUN_SBE_H
#define SBE_WH_CU_RUN_SBE_H

#include "runSBE.h"

#include "CuWhFourierTransform.h"
#include "CuDiagonalizator.h"
#include "CuPropagator.h"
#include "CuObserver.h"

#include <boost/numeric/odeint/external/thrust/thrust.hpp>

int runSBE_GPU(Parameter_t &param);

#endif
