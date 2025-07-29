#include "PureGauss1D.h"

PureGauss1D::PureGauss1D(double E0, double FWHM, double tStart_FWHM, double tEnd_FWHM) :
    FWHM(FWHM),
    Pulse1D(E0, tStart_FWHM * FWHM, tEnd_FWHM * FWHM)
{
    beta =  4 * log(2) / (FWHM * FWHM);
    Aprefactor = E0 * sqrt( M_PI / beta) / 2;
    erfStartTime = Aprefactor * erf ( getStartTime() );
}

PureGauss1D::~PureGauss1D()
{

}

double PureGauss1D::getA0() const
{
    return -getA(endTime);
}

double PureGauss1D::getRisingFWHM() const
{
    return FWHM;
}

double PureGauss1D::getE(double t) const
{
    return E0 * exp(-beta * t * t);
}

double PureGauss1D::getA(double t) const
{
    return - (Aprefactor * erf( t * sqrt(beta) ) - erfStartTime);
}

const char * PureGauss1D::getHumanReadableType() const {
    return "pure gaussian";
}

void PureGauss1D::print(FILE *f) const
{
    int indent = sbewhConfig::pulseParameterPadLength;
    Pulse1D::print(f);
    fprintf(f, "%-*s%.1lf\t (%.2lf fs)\n", indent, "FWHM", FWHM,
               atomicUnits::to_fs(FWHM));
}
