#include "PropagatorBase.h"

PropagatorBase::PropagatorBase(const StateContext &sc) : sc(sc)
{

}

void PropagatorBase::timeAndPrint(const std::chrono::steady_clock::time_point &startTime)
{
    std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
    allTime_us += (unsigned)std::chrono::duration_cast<std::chrono::microseconds>(endTime - lastCallTime).count();
    derivTime_us += (unsigned)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
    callCount++;
    if ( callCount % printEveryCalls == 0){
        callCount = 0;
        Logger::info("Averaged time for:\n\t\tComplete stepping: %u us\n\t\tPure dsdt: %u us\n",
                              allTime_us / printEveryCalls, derivTime_us / printEveryCalls);
        allTime_us = 0;
        derivTime_us = 0;
    }
    lastCallTime = std::chrono::steady_clock::now();
}

void PropagatorBase::reset()
{
    lastCallTime = std::chrono::steady_clock::now();
    allTime_us = 0;
    derivTime_us = 0;
    callCount = 0;
}
