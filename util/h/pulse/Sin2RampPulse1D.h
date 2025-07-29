#ifndef SBE_WH_SIN2_RAMP_PULSE_H
#define SBE_WH_SIN2_RAMP_PULSE_H


#include "math.h"
#include "Pulse1D.h"

struct Sin2RampPulse1Draw_t {
    double beta;
    double omega;
    double risingCycles;
    int cyclesOn;
    Pulse1Draw_t pulse1D;
};
/*
 * Defines oszillating pulse with envelope given by sin^2(FWHM) - constant - sin^2(FWHM)  in 1D
 */
class Sin2RampPulse1D : public Pulse1D {
private:
    double omega;
    double risingCycles;
    int cyclesOn;
    double beta;
public:
    Sin2RampPulse1D(double E0, double omega, double risingCycles, int cyclesOn);
    Sin2RampPulse1D(const Sin2RampPulse1Draw_t & r) : Pulse1D(r.pulse1D), omega(r.omega), risingCycles(r.risingCycles),
                                                                cyclesOn(r.cyclesOn), beta(r.beta) {}
    ~Sin2RampPulse1D();
    Sin2RampPulse1Draw_t getRaw() const { return Sin2RampPulse1Draw_t { .beta = beta, .omega = omega, .risingCycles = risingCycles,
                                                                                  .cyclesOn = cyclesOn, .pulse1D = Pulse1D::getBaseRaw()}; }
    double getA0() const override;
    double getRisingFWHM() const override;
    double getE(double t) const override;
    double getA(double t) const override;
    double getOmega() const;
    double getCyclesOn() const;

    const char * getHumanReadableType() const override;
    Pulse1Dtype getType() const override { return Sin2RampPulse; }
    void print(FILE *f) const override;
};

#endif
