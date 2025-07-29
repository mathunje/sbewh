#include "Sin2Pulse1D.h"

Sin2Pulse1D::Sin2Pulse1D(double E0, double omega, double cep, double FWHM, double tStart_FWHM, double tEnd_FWHM) :
    FWHM(FWHM),
    omega(omega),
    cep(cep),
    Pulse1D(E0, tStart_FWHM * FWHM, tEnd_FWHM * FWHM)
{
    beta = M_PI / 2 / FWHM;
}

Sin2Pulse1D::~Sin2Pulse1D()
{

}

double Sin2Pulse1D::getA0() const
{
    return E0 / omega;
}

double Sin2Pulse1D::getRisingFWHM() const
{
    return FWHM;
}

double Sin2Pulse1D::getE(double t) const
{
    if ( t < -FWHM || t > FWHM)
        return 0;
    double cp = omega * t + cep;
    double envCos = cos(beta * t);
    return E0 *(cos(cp)*envCos*envCos - 2 * beta / omega  * sin(cp)*envCos*sin(beta*t) );
}

double Sin2Pulse1D::getA(double t) const
{
    if ( t < -FWHM || t > FWHM)
        return 0;
    double envCos = cos(beta * t);
    return - E0 / omega  * envCos * envCos * sin(omega * t + cep);
}

const char * Sin2Pulse1D::getHumanReadableType() const {
    return "sin2 pulse";
}

double Sin2Pulse1D::getOmega() const
{
    return omega;
}

double Sin2Pulse1D::getCep() const
{
    return cep;
}

void Sin2Pulse1D::print(FILE *f) const
{
    int indent = sbewhConfig::pulseParameterPadLength;
    Pulse1D::print(f);
    fprintf(f, "%-*s%.1lf\t (%.2lf fs)\n", indent, "FWHM", FWHM,
               atomicUnits::to_fs(FWHM));
    fprintf(f, "%-*s%.5lf\n", indent, "cep", cep);
    fprintf(f, "%-*s%.5lf\t (%.2lf nm)\n", indent, "lambda", omega,
                atomicUnits::to_nm( atomicUnits::speedOfLight * 2 * M_PI / omega));
}
