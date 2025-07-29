#include "io/GeneralParameter.h"

GeneralParameterValidator::GeneralParameterValidator(ParameterIO &paramIO)
{
    std::map<unsigned, std::string> description;
    rmMap["SBE"] = RunMode::SBE;
    description[(unsigned)(RunMode::SBE)] = "propagation of density matrix";
    rmMap["WHtransform"] = RunMode::WHtransform;
    description[(unsigned)(RunMode::WHtransform)] = "writes matrix elements (Wannier basis) on k-grid into a file";
    rmMap["FermiLevel"] = RunMode::FermiLevel;
    description[(unsigned)(RunMode::FermiLevel)] = "calculates Fermi level with more details";
    rmMap["Pulse"] = RunMode::Pulse;
    description[(unsigned)(RunMode::Pulse)] = "outputs pulse with its subpulse contributions into a file";

    std::string mergedDescription;
    std::string allowedModes;
    for(auto it = rmMap.begin(); it != rmMap.end(); ++it){
        if (it == rmMap.begin()){
            allowedModes += it->first;
        } else {
            allowedModes += "," + it->first;
        }
        const std::string &subDes = description[(unsigned)(it->second)];
        mergedDescription += "  " + it->first + ":\n"
                           + "    " + subDes + "\n";
    }
    paramIO.startGroup("General", "comments on the input format:\n"
                                  "- vectors of numbers may be separated by space or colon\n"
                                  "- lines maybe split by using a backslash");
    paramIO.registerParam(new ParameterString("runMode", "defines run mode (" + allowedModes + ")\n"
                                                          + mergedDescription, runModeStr, "Dry"));
    paramIO.registerParam(new ParameterBool("dryRun", "only parses parameter", dryRun, false));
    ParameterNumeric<unsigned> * pDim = new ParameterNumeric<unsigned>("dimensionality",
                                                    "defines, whether 1D, 2D or 3D calculation is performed", dim, 3);
    pDim->setRange(1, 3);
    paramIO.registerParam(pDim);
    paramIO.registerParam(new ParameterBool("verbose", "if set, output is more verbose", verbose, false));
    paramIO.registerParam(new ParameterBool("veryVerbose", "if set, output is very verbose", veryVerbose, false));
    paramIO.registerParam(new ParameterBool("timeIt", "prints some additional timing information", timeIt, false));
    paramIO.registerParam(new ParameterBool("finishMultiRuns", "starts unfinished (and not already running) multi runs of last/specified run directory\n"
                                                               "Intended to be used only as command line argument\n"
                                                               "All other parameters are ignored\n", finishMultiRuns, false));
    paramIO.finishGroup();
}

bool GeneralParameterValidator::update(GeneralParameter_t *p)
{
    auto it = rmMap.find(runModeStr);
    if ( it == rmMap.end()){
        Logger::error("Unrecognized run mode '%s'\n", runModeStr.c_str());
        return false;
    } else {
        p->runMode = it->second;
    }
    p->dryRun = dryRun;
    p->dim = dim;
    p->verbose = verbose || veryVerbose;
    p->veryVerbose = veryVerbose;
    p->timeIt = timeIt;
    p->finishMultiRuns = finishMultiRuns;
    return true;
}

void mpi_bcast_GeneralParameter(GeneralParameter_t &p, const MpiParameter_t &mp)
{
    mpi_bcast_trivial(p.dryRun, mp);
    mpi_bcast_trivial(p.dim, mp);
    mpi_bcast_trivial(p.verbose, mp);
    mpi_bcast_trivial(p.veryVerbose, mp);
    mpi_bcast_trivial(p.timeIt, mp);
    mpi_bcast_trivial(p.runMode, mp);
    mpi_bcast_trivial(p.finishMultiRuns, mp);
}
