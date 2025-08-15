#ifndef SBE_WH_GAUSSIAN_PULSE_1D_H
#define SBE_WH_GAUSSIAN_PULSE_1D_H


#include "math.h"
#include "Pulse1D.h"


struct GaussianPulse1Draw_t {
    double beta;
    double FWHM;
    double omega;
    double cep;
    Pulse1Draw_t pulse1D;
};

/*
 * Defines oszillating pulse with gaussian envelope in 1D
 */
class GaussianPulse1D : public Pulse1D {
private:
    double FWHM;
    double omega;
    double cep;
    double beta;
public:
    GaussianPulse1D(double E0, double omega, double cep, double FWHM, double tStart_FWHM=-3, double tEnd_FWHM = 3);
    GaussianPulse1D(const GaussianPulse1Draw_t & r) : Pulse1D(r.pulse1D), FWHM(r.FWHM), omega(r.omega), cep(r.cep), beta(r.beta) {}
    ~GaussianPulse1D();
    GaussianPulse1Draw_t getRaw() const { return GaussianPulse1Draw_t { .beta = beta, .FWHM = FWHM, .omega = omega, .cep = cep, .pulse1D = Pulse1D::getBaseRaw()}; }
    double getA0() const override;
    double getRisingFWHM() const override;
    double getE(double t) const override;
    double getA(double t) const override;
    double getOmega() const;
    double getCep() const;

    const char * getHumanReadableType() const override;
    Pulse1Dtype getType() const override { return Gaussian; }
    void print(FILE *f) const override;
};

#endif
