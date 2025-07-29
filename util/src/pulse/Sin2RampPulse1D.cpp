#include "Sin2RampPulse1D.h"

Sin2RampPulse1D::Sin2RampPulse1D(double E0, double omega, double risingCycles, int cyclesOn) :
    omega(omega),
    risingCycles(risingCycles),
    cyclesOn(cyclesOn),
    Pulse1D(E0, - 2 * M_PI / omega * (risingCycles + 0.5*cyclesOn), 2 * M_PI / omega * (risingCycles + 0.5 * cyclesOn) )
{
    beta = omega / 4 / risingCycles;
}

Sin2RampPulse1D::~Sin2RampPulse1D()
{

}

double Sin2RampPulse1D::getA0() const
{
    return E0 / omega;
}

double Sin2RampPulse1D::getRisingFWHM() const
{
    return M_PI / (2 * beta);
}


double Sin2RampPulse1D::getE(double t) const
{
    double risingTime = 2 * M_PI / omega * risingCycles;
    double onTime = M_PI / omega * cyclesOn;
    double halfTime = risingTime + M_PI / omega * cyclesOn;
    if ( t < -halfTime )
        return 0;
    if ( t < -onTime ){
        t += onTime;
        double cp = omega * t;
        double envCos = cos(beta * t);
        return E0 *(cos(cp)*envCos*envCos - 2 * beta / omega  * sin(cp)*envCos*sin(beta*t) );
    }
    if ( t < onTime ){
        return ((cyclesOn % 2 == 0) ? 1 : -1) * E0 * cos(omega * t);
    }
    if ( t < halfTime){
        t -= onTime;
        double cp = omega * t;
        double envCos = cos(beta * t);
        return E0 *(cos(cp)*envCos*envCos - 2 * beta / omega  * sin(cp)*envCos*sin(beta*t) );
    }
    return 0;
}

double Sin2RampPulse1D::getA(double t) const
{

    double risingTime = 2 * M_PI / omega * risingCycles;
    double onTime = M_PI / omega * cyclesOn;
    double halfTime = risingTime + M_PI / omega * cyclesOn;
    if ( t < -halfTime )
        return 0;
    if ( t < -onTime ){
        t += onTime;
        double envCos = cos(beta * t);
        return - E0 / omega  * envCos * envCos * sin(omega * t);
    }
    if ( t < onTime ){
        return ((cyclesOn % 2 == 0) ? 1 : -1) * -E0 / omega * sin(omega * t);
    }
    if ( t < halfTime){
        t -= onTime;
        double envCos = cos(beta * t);
        return - E0 / omega  * envCos * envCos * sin(omega * t);
    }
    return 0;
}

const char * Sin2RampPulse1D::getHumanReadableType() const {
    return "sin2Ramp pulse";
}

double Sin2RampPulse1D::getOmega() const
{
    return omega;
}

double Sin2RampPulse1D::getCyclesOn() const
{
    return cyclesOn;
}

void Sin2RampPulse1D::print(FILE *f) const
{
    int indent = sbewhConfig::pulseParameterPadLength;
    Pulse1D::print(f);
    fprintf(f, "%-*s%.1lf\n", indent, "rising cycles", risingCycles);
    fprintf(f, "%-*s%d\n", indent, "cycles on", cyclesOn);
    fprintf(f, "%-*s%.5lf\t (%.2lf nm)\n", indent, "lambda", omega,
                atomicUnits::to_nm( atomicUnits::speedOfLight * 2 * M_PI / omega));
}
