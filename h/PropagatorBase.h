#ifndef SBE_WH_PROPAGATOR_BASE_H
#define SBE_WH_PROPAGATOR_BASE_H

#include "Parameter.h"
#include "StateContext.h"

#include <cmath>
#include <chrono>

class PropagatorBase{
protected:
    const StateContext &sc;
    std::chrono::steady_clock::time_point lastCallTime;
    unsigned allTime_us;
    unsigned derivTime_us;
    unsigned callCount;
    const unsigned printEveryCalls = 1000;
    void timeAndPrint(const std::chrono::steady_clock::time_point &startTime);
public:
    PropagatorBase(const StateContext &sc);
    void reset();
};

#endif
