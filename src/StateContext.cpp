#include "StateContext.h"

GeomVector3d kFracTo_au(const GeomVector3d &kFrac, const std::array<GeomVector3d, 3> &Rl)
{
    double preFact = 2 * std::numbers::pi / ( Rl[0].dot(Rl[1].cross(Rl[2])));
    return preFact * ( kFrac.at(0) * Rl[1].cross(Rl[2]) +
                       kFrac.at(1) * Rl[2].cross(Rl[0]) +
                       kFrac.at(2) * Rl[0].cross(Rl[1]) );
}


StateContext::StateContext(const Parameter_t &param, unsigned NrProcStart, unsigned NrProcCount) :
    param(param),
    kGlobalShift_au(kFracTo_au(param.prop.kGlobalShift, param.tb.latticeVectors)),
    NrProcStart(NrProcStart),
    NrProcCount(NrProcCount),
    numWann(param.tb.numWann),
    dim(param.general.dim),
    si(numWann)
{

}
