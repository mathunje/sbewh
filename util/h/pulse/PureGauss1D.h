#ifndef SBE_WH_PULSE_PURE_GAUSS_1D_H
#define SBE_WH_PULSE_PURE_GAUSS_1D_H

#include "math.h"
#include "Pulse1D.h"


struct PureGauss1Draw_t {
    double beta;
    double FWHM;
    double erfStartTime;
    double Aprefactor;
    Pulse1Draw_t pulse1D;
};

/*
 * Defines a non-oszillating gaussian beam in 1D
 * It is of course unphysical, but may be used to calculate absorption spectra
 */
class PureGauss1D : public Pulse1D {
private:
    double FWHM;
    double beta;
    double erfStartTime;
    double Aprefactor;
public:
    PureGauss1D(double E0, double FHWM, double tStart_FHWM = -3, double tEnd_FHWM = 3);
    PureGauss1D(const PureGauss1Draw_t & r) : Pulse1D(r.pulse1D), FWHM(r.FWHM), beta(r.beta), erfStartTime(r.erfStartTime), Aprefactor(r.Aprefactor){}
    ~PureGauss1D();
    PureGauss1Draw_t getRaw() const { return PureGauss1Draw_t { .beta = beta, .FWHM = FWHM, .erfStartTime = erfStartTime, .Aprefactor = Aprefactor,
                                                                              .pulse1D = Pulse1D::getBaseRaw()}; }
    double getA0() const override;
    double getRisingFWHM() const override;
    double getE(double t) const override;
    double getA(double t) const override;
    const char * getHumanReadableType() const override;
    Pulse1Dtype getType() const override { return PureGauss; }
    void print(FILE *f) const override;
};

#endif
