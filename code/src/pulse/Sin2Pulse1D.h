#ifndef SBE_WH_SIN2_PULSE_1D_H
#define SBE_WH_SIN2_PULSE_1D_H

#include "math.h"
#include "Pulse1D.h"

struct Sin2Pulse1Draw_t {
    double beta;
    double FWHM;
    double omega;
    double cep;
    Pulse1Draw_t pulse1D;
};
/*
 * Defines oszillating pulse with sin^2 envelope in 1D
 */
class Sin2Pulse1D : public Pulse1D {
private:
    double FWHM;
    double omega;
    double cep;
    double beta;
public:
    Sin2Pulse1D(double E0, double omega, double cep, double FWHM, double tStart_FWHM=-1, double tEnd_FWHM = 1);
    Sin2Pulse1D(const Sin2Pulse1Draw_t & r) : Pulse1D(r.pulse1D), FWHM(r.FWHM), omega(r.omega), cep(r.cep), beta(r.beta) {}
    ~Sin2Pulse1D();
    Sin2Pulse1Draw_t getRaw() const { return Sin2Pulse1Draw_t { .beta = beta, .FWHM = FWHM, .omega = omega, .cep = cep, .pulse1D = Pulse1D::getBaseRaw()}; }
    double getA0() const override;
    double getRisingFWHM() const override;
    double getE(double t) const override;
    double getA(double t) const override;
    double getOmega() const;
    double getCep() const;

    const char * getHumanReadableType() const override;
    Pulse1Dtype getType() const override { return Sin2Pulse; }
    void print(FILE *f) const override;
};

#endif
