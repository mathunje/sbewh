#ifndef SBE_WH_STATE_CONTEXT_H
#define SBE_WH_STATE_CONTEXT_H

#include "Parameter.h"
#include "GeomVector.hpp"

struct StateIndexMapper{
private:
    const unsigned padSI = 8;
    unsigned cid = 0;
    unsigned nextId() { return cid++; }
    inline unsigned roundUp(unsigned v) { return padSI * ((v+padSI-1)/padSI); }
public:
    const unsigned paddedMatSize;
    const unsigned A[3];
    const unsigned specialIndexCount;

    StateIndexMapper(unsigned numWann) :
        paddedMatSize( roundUp(numWann*numWann) ),
        A {nextId(), nextId(), nextId()},
        specialIndexCount( roundUp(nextId()) )
        {}

    inline unsigned rhoK(unsigned id) const { return specialIndexCount + paddedMatSize * id; }  // starting index of matrix rho_k
};

struct StateContext {
    StateContext(const Parameter_t &param, unsigned NrProcStart=0, unsigned NrProcCount=0);
    const Parameter_t &param;
    const GeomVector3d kGlobalShift_au;
    unsigned NrProcStart;
    const unsigned NrProcCount;
    const unsigned numWann;
    const unsigned dim;
    const StateIndexMapper si;
    GeomVector3d getE(double t) const { return param.pulse.pulse.getE(t); }
};


#endif
