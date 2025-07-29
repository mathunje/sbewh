#include "Pulse1D.h"

Pulse1D::Pulse1D(double E0, double startTime, double endTime) :
        E0(E0), startTime(startTime), endTime(endTime){
}

Pulse1D::~Pulse1D()
{

}

double Pulse1D::getStartTime() const
{
    return startTime;
}

double Pulse1D::getEndTime() const
{
    return endTime;
}

double Pulse1D::getE0() const
{
    return E0;
}

void Pulse1D::print(FILE *f) const
{
    int indent = sbewhConfig::pulseParameterPadLength;
    fprintf(f, "Pulse of type %s\n", getHumanReadableType());
    fprintf(f, "%-*s%.4lf\t (%.2lf V/nm)\t (-> %.3e W/cm^2)\n", indent, "E0", E0,
                   atomicUnits::to_V_nm(E0), atomicUnits::to_W_cm2(E0*E0));
    fprintf(f, "%-*s%.1lf\t (%.2lf fs)\n", indent, "startTime", startTime,
                   atomicUnits::to_fs(startTime));
    fprintf(f, "%-*s%.1lf\t (%.2lf fs)\n", indent, "endTime", endTime,
                   atomicUnits::to_fs(endTime));

}

double Pulse1D::numericalIntegrationE() const{
    return numericalIntegrationE(startTime, endTime);
}

double Pulse1D::numericalIntegrationE(double startT, double endT) const
{
    auto E = [&] (double t) { return getE(t); };
    return boost::math::quadrature::trapezoidal(E, startT, endT);
}


double Pulse1D::numericalIntegrationErrorE() const
{
    auto error = [&] (double t) { double diff = numericalIntegrationE(startTime, t)  + getA(t);
                                  return diff * diff; };
    return boost::math::quadrature::trapezoidal(error, startTime, endTime);
}
