#include "Pulse3D.h"

Pulse3D::Pulse3D()
{

}

void Pulse3D::addPulse(const GeomVector3d &pol, std::shared_ptr<Pulse1D> pulse, double centeredAt)
{
    GeomVector3d npol = pol;
    npol.normalize();
    pulses.push_back({npol, pulse, centeredAt});
}

GeomVector3d Pulse3D::getE(double t) const
{
    GeomVector3d E;
    for(const auto &p : pulses)
        E += p.pulse1D->getE(t - p.center) * p.pol;
    return E;
}

GeomVector3d Pulse3D::getA(double t) const
{
    GeomVector3d A;
    for(const auto &p : pulses)
        A += p.pulse1D->getA(t - p.center) * p.pol;
    return A;
}

GeomVector3d Pulse3D::getSinglePulseE(unsigned i, double t) const
{
    return pulses[i].pulse1D->getE(t - pulses[i].center) * pulses[i].pol;
}

GeomVector3d Pulse3D::getSinglePulseA(unsigned i, double t) const
{
    return pulses[i].pulse1D->getA(t - pulses[i].center) * pulses[i].pol;
}


double Pulse3D::getStartTime() const
{
    double startTime = std::numeric_limits<double>::infinity();
    for(const auto &p : pulses){
        double t = p.pulse1D->getStartTime() + p.center;
        if ( t < startTime )
            startTime = t;
    }
    return startTime;
}

double Pulse3D::getEndTime() const
{
    double endTime = -std::numeric_limits<double>::infinity();
    for(const auto &p : pulses){
        double t = p.pulse1D->getEndTime() + p.center;
        if ( t > endTime )
            endTime = t;
    }
    return endTime;
}

void Pulse3D::clear()
{
    pulses.clear();
}


void Pulse3D::print(FILE *f) const
{
    int indent = sbewhConfig::pulseParameterPadLength;
    fprintf(f, "3D pulse composed of:\n");
    for(const auto &p: pulses){
        p.pulse1D->print(f);
        fprintf(f, "%-*s%.1lf\t (%.2lf fs)\n", indent, "center", p.center,
                   atomicUnits::to_fs(p.center));
        fprintf(f, "%-*s(%lf, %lf, %lf)\n", indent, "polarization", p.pol.at(0), p.pol.at(1), p.pol.at(2) );
        fprintf(f, "\n");
    }
}
