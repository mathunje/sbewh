#include "io/OutputParameter.h"

OutputParameterValidator::OutputParameterValidator(ParameterIO &paramIO){
    paramIO.startGroup("Output", "defined if/how output files should be written");
    paramIO.registerParam(new ParameterBool("disableRunSaving", "disables saving of run results", disableSaving, false));
    paramIO.registerParam(new ParameterBool("saveInput", "if set: input files are saved at program start", saveInput, true));
    auto * pswi = new ParameterBool("saveWannierInput", "if set: wannier files are saved at program start"
                                                        "(overrides saveInput for Wannier files)", saveWannierInput, true);
    saveWannierInputSetPtr = pswi->getSetPtr();
    paramIO.registerParam(pswi);
    paramIO.registerParam(new ParameterString("outputDirectory", "output directory where runs are stored", baseDir, "runs"));
    paramIO.registerParam(new ParameterBool("saveAsNpz", "uses purely npz files for output (otherwise some *.txt files are generted)", saveAsNpz, true));
    paramIO.registerParam(new ParameterString("run", "if set defines run directory otherwise the run directory is choosen automatically.", run, ""));
    paramIO.finishGroup();
}

bool OutputParameterValidator::update(OutputParameter_t *p, bool selectLastRun){
    p->disableSaving = disableSaving;
    p->saveInput = saveInput && ! disableSaving;
    p->saveAsNpz = saveAsNpz;
    if ( *saveWannierInputSetPtr ){
        p->saveWannierInput = saveWannierInput && ! disableSaving;
    } else {
        p->saveWannierInput = p->saveInput;
    }
    if ( p->saveWannierInput && ! p->saveInput){
        Logger::warn("Saving input (even if not requested) because Wannier input was requested\n");
        p->saveInput = true;
    }
    p->baseDir = baseDir;
    if ( disableSaving )
        return true;
    namespace fs = std::filesystem;
    auto bd = fs::path(baseDir);
    if ( run.compare("") ){
        if ( fs::is_directory(fs::status(bd / run) ) ){
            if ( ! selectLastRun ){
                Logger::error("Run '%s' is already existent in '%s'\n", run.c_str(), baseDir.c_str() );
                return false;
            }
        } else {
            if ( selectLastRun ){
                Logger::error("Run '%s' was not setup yet in '%s'\n", run.c_str(), baseDir.c_str());
                return false;
            }
        }
    } else {
        if ( fs::is_directory(fs::status(bd)) ){
            char runName[8];
            int runNumber = 0;
            do {
                sprintf(runName, "%03d", ++runNumber);
            } while ( fs::exists(fs::status(bd / fs::path(runName) )));
            if ( selectLastRun ){
                runNumber--;
                if (runNumber < 1){
                    Logger::error("output directory '%s' contains no automatically chosen sub run\n", baseDir.c_str());
                    return false;
                }
            }
            sprintf(runName, "%03d", runNumber);
            run = std::string(runName);
        } else {
            if ( selectLastRun ){
                Logger::error("output directory '%s' not yet existent\n", baseDir.c_str());
                return false;
            }
            run = "001";
        }
    }
    p->run = run;
    p->saveDir = bd / run;
    return true;
}

void mpi_bcast_OutputParameter(OutputParameter_t &p, const MpiParameter_t &mp)
{
    mpi_bcast_trivial(p.disableSaving, mp);
    mpi_bcast_trivial(p.saveInput, mp);
    mpi_bcast_trivial(p.saveWannierInput, mp);
    mpi_bcast_trivial(p.saveAsNpz, mp);
    mpi_bcast(p.baseDir, mp);
    mpi_bcast(p.saveDir, mp);
    mpi_bcast(p.run, mp);
}
