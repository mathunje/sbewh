#ifndef SBE_WH_PULSE_3D_H
#define SBE_WH_PULSE_3D_H

#include "sbewhConfig.h"

#include <limits>
#include <vector>
#include <memory>

#include <stdio.h>

#include "Pulse1D.h"
#include "util/GeomVector.hpp"


struct SinglePulse3D_t {
    GeomVector3d pol;
    std::shared_ptr<Pulse1D> pulse1D;
    double center;
};

class Pulse3D {
    std::vector<SinglePulse3D_t> pulses;
public:
    Pulse3D();
    void addPulse(const GeomVector3d &pol, std::shared_ptr<Pulse1D> pulse, double centeredAt=0);
    void addPulse(const SinglePulse3D_t s) { pulses.push_back(s); }
    const SinglePulse3D_t & getSinglePulse(unsigned id) const { return pulses[id]; }
    unsigned getSinglePulseCount() const { return pulses.size(); }
    void clear();
    double getStartTime() const;
    double getEndTime() const;
    GeomVector3d getE(double t) const;
    GeomVector3d getA(double t) const;
    GeomVector3d getSinglePulseE(unsigned i, double t) const;
    GeomVector3d getSinglePulseA(unsigned i, double t) const;
    void print(FILE *f) const;
};

#endif
