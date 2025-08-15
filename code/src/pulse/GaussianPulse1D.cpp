#include "GaussianPulse1D.h"

GaussianPulse1D::GaussianPulse1D(double E0, double omega, double cep, double FWHM, double tStart_FWHM, double tEnd_FWHM) :
    FWHM(FWHM),
    omega(omega),
    cep(cep),
    Pulse1D(E0, tStart_FWHM * FWHM, tEnd_FWHM * FWHM)
{
    beta = 4 * log(2) / (FWHM * FWHM);
}

GaussianPulse1D::~GaussianPulse1D()
{

}

double GaussianPulse1D::getA0() const
{
    return E0 / omega;
}

double GaussianPulse1D::getRisingFWHM() const
{
    return FWHM;
}

double GaussianPulse1D::getE(double t) const
{
    double cphase = omega * t + cep;
    return E0 * ( - 2 * beta/omega  * t * sin (cphase) + cos ( cphase )) * exp(-beta * t * t);
}

double GaussianPulse1D::getA(double t) const
{
    return - E0 / omega * exp(-beta * t * t) * sin ( omega * t + cep);
}

const char * GaussianPulse1D::getHumanReadableType() const {
    return "gaussian pulse";
}

double GaussianPulse1D::getOmega() const
{
    return omega;
}

double GaussianPulse1D::getCep() const
{
    return cep;
}

void GaussianPulse1D::print(FILE *f) const
{
    int indent = sbewhConfig::pulseParameterPadLength;
    Pulse1D::print(f);
    fprintf(f, "%-*s%.1lf\t (%.2lf fs)\n", indent, "FWHM", FWHM,
               atomicUnits::to_fs(FWHM));
    fprintf(f, "%-*s%.5lf\n", indent, "cep", cep);
    fprintf(f, "%-*s%.5lf\t (%.2lf nm)\n", indent, "lambda", omega,
                atomicUnits::to_nm( atomicUnits::speedOfLight * 2 * M_PI / omega));
}
