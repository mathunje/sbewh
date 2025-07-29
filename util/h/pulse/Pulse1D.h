#ifndef SBE_WH_PULSE_1D_H
#define SBE_WH_PULSE_1D_H

#include "sbewhConfig.h"

#include <stdio.h>
#include <vector>
#include<boost/math/quadrature/trapezoidal.hpp>
#include "unitConversion.h"


enum Pulse1Dtype { Gaussian, PureGauss, Sin2Pulse, Sin2RampPulse};

struct Pulse1Draw_t {
    double startTime;
    double endTime;
    double E0;
};

/*
 * Defines Baseclass for all 1D pulse with basic support for
    - E field integration
    - printing
    - plotting
 */
class Pulse1D {
protected:
    double startTime;
    double endTime;
    double E0;
public:
    Pulse1D(double E0, double startTime, double endTime);
    Pulse1D(const Pulse1Draw_t &r) : startTime(r.startTime), endTime(r.endTime), E0(r.E0)  {}
    virtual ~Pulse1D();
    virtual Pulse1Draw_t getBaseRaw() const { return Pulse1Draw_t { .startTime = startTime, .endTime = endTime, .E0 = E0}; }
    virtual double getStartTime() const final;
    virtual double getEndTime() const final;
    virtual double getE0() const final;
    virtual double getA0() const = 0; // return resonable estimate for maximum modulus of A.
    virtual double getRisingFWHM() const = 0;
    virtual double getE(double t) const = 0;
    virtual double getA(double t) const = 0;

    virtual const char * getHumanReadableType() const = 0;
    virtual Pulse1Dtype getType() const = 0;
    virtual void print(FILE *f) const;

    double numericalIntegrationE() const; // integral E(t)dt from startTime to endTime
    double numericalIntegrationE(double startT, double endT) const; // integral E(t)dt

    double numericalIntegrationErrorE() const; // integral | integral E(t')dt' from startTime to t + A(t) |^2 dt
};

#endif
