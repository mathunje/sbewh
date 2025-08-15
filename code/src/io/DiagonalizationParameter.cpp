#include "io/DiagonalizationParameter.h"


DiagonalizationParameterValidator::DiagonalizationParameterValidator(ParameterIO &paramIO)
{
    dmMap["jacobi"] = DiagMode::jacobi;
#ifdef SBE_WH_LAPACK
    dmMap["zheev"] = DiagMode::zheev;
    dmMap["zheevd"] = DiagMode::zheevd;
    for(auto it = dmMap.begin(); it != dmMap.end(); ++it){
        if (it == dmMap.begin()){
            allowedModes += it->first;
        } else {
            allowedModes += "," + it->first;
        }
    }
    const char * defaultDiagMode = "zheev";
#else
    const char * defaultDiagMode = "jacobi";
#endif
    paramIO.startGroup("Diagonalization", "Select lapack routine or threshold for jacobi diagonalization");
    paramIO.registerParam(new ParameterString("mode", "defines diagonalization routine (" + allowedModes + ")", diagModeStr, defaultDiagMode));
    ParameterNumeric<unsigned> *pMSC = new ParameterNumeric<unsigned>("maxSweepCount",
                                                                      "maximum number of jacobi sweeps",
                                                                      maxSweepCount, 5);
    pMSC->setRange(1, 50);
    paramIO.registerParam(pMSC);
    ParameterNumeric<double> *pMRE = new ParameterNumeric<double>("maxRelError",
                                                                  "maximum allowed relative error of jacobi diagonalization",
                                                                   maxRelError, 1e-14);
    pMRE->setRange(1e-15, 1e-2);
    paramIO.registerParam(pMRE);
    paramIO.finishGroup();
}

bool DiagonalizationParameterValidator::update(DiagonalizationParameter_t *p)
{
    auto it = dmMap.find(diagModeStr);
    if ( it == dmMap.end()){
        Logger::error("Unrecognized diagonalization mode '%s'\n", diagModeStr.c_str());
        return false;
    }
    p->mode = it->second;
    p->maxSweepCount = maxSweepCount;
    p->maxRelError = maxRelError;
    return true;
}

void mpi_bcast_DiagonalizationParameter(DiagonalizationParameter_t &p, const MpiParameter_t &mp)
{
    mpi_bcast_trivial(p.mode, mp);
    mpi_bcast_trivial(p.maxSweepCount, mp);
    mpi_bcast_trivial(p.maxRelError, mp);
}
