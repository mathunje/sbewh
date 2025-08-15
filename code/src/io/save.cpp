#include "io/save.h"


void saveInput(const Parameter_t &param)
{
    // save parameter input file
    namespace fs = std::filesystem;
    fs::path p(param.out.saveDir);
    for(const auto & fname : param.cmdInp.fnames)
        fs::copy(fname, p / fs::path(fname).filename() );
    FILE * fArgs = fopen(( p / fs::path("args.txt")).c_str(), "w");
    fprintf(fArgs, "Git version: %s\n", kGitHash);
    fprintf(fArgs, "Command line arguments (one per line):\n");
    for(unsigned i=0; i<param.cmdInp.args.size(); i++)
        fprintf(fArgs, "%s\n", param.cmdInp.args[i].c_str());
    fclose(fArgs);
    // save input files
    if ( param.out.saveWannierInput){
        for (auto fn : { param.tb.wannierTbFname, param.tb.wannierWsvecFname})
            if ( fn.size() && fs::exists(fs::status(fn)) )
                fs::copy(fn, p / fs::path(fn).filename() );
    }
    for (auto fn : { param.pulse.multiRunFname})
        if ( fn.size() && fs::exists(fs::status(fn)) )
            fs::copy(fn, p / fs::path(fn).filename() );
}


inline bool saveRegionSample(const Parameter_t & param, const std::filesystem::path &p, std::string fname, const ExpectationValues_t &ev,
                             unsigned rid, unsigned sampleIndex)
{
    unsigned numWann = param.tb.numWann;
    const RegionExpectationValues_t & rexp = ev.region[rid];
    unsigned dim = rexp.dim[0] * rexp.dim[1] * rexp.dim[2];
    FILE * fev = fopen( ( p / fname).c_str() , "w");
    if ( ! fev ){
        Logger::error("Could not open '%s'\n", fname.c_str() );
        return false;
    }
    fprintf(fev, "# all in atomic units\n");
    fprintf(fev, "# t\tweight\tJx\tJy\tJz\t");
    for(unsigned bi=0; bi<numWann; bi++)
        fprintf(fev, "\toccW%u", bi);
    for(unsigned bi=0; bi<numWann; bi++)
        fprintf(fev, "\toccH%u", bi);
    fprintf(fev, "\n");
    for(unsigned ti=0; ti<ev.t.size(); ti++){
        unsigned ts = ti * dim + sampleIndex;
        fprintf(fev, "%lg\t%lg", ev.t[ti], rexp.kPointWeightSum[ts]);
        for(unsigned dir=0; dir<3; dir++)
            fprintf(fev, "\t%lg", rexp.j[ts][dir]);
        for(unsigned i=ts*numWann; i<(ts+1)*numWann; i++)
            fprintf(fev, "\t%.10lg", rexp.occupationWan[i]);
        for(unsigned i=ts*numWann; i<(ts+1)*numWann; i++)
            fprintf(fev, "\t%.10lg", rexp.occupationHam[i]);
        fprintf(fev, "\n");
    }
    if ( fclose(fev) ) {
        Logger::error("Expectation value file '%s' not written properly\n", fname.c_str());
        return false;
    }

    if ( param.ksr.storeDensityMatrix ){
        std::string rhoName = p / ( "rho_" + fname);
        FILE * fev = fopen( rhoName.c_str(), "w");
        if ( ! fev ){
            Logger::error("Could not open '%s'\n", rhoName.c_str() );
            return false;
        }
        fprintf(fev, "# interpolated density matrix, time in atomic units\n");
        fprintf(fev, "# t\tweight");
        for(unsigned i=0; i<numWann; i++)
            for(unsigned j=0; j<numWann; j++)
                fprintf(fev, "\tRe_%u_%u\tIm_%u_%u", i, j, i, j);
        fprintf(fev, "\n");
        for(unsigned tdi=0; tdi<ev.tDensIndices.size(); tdi++){
            double ti = ev.tDensIndices[tdi];
            unsigned tds = tdi * dim + sampleIndex;
            fprintf(fev, "%lg\t%lg", ev.t[ti], rexp.kPointWeightSum[tds]);
            for(unsigned i=tds*numWann*numWann; i<(tds+1)*numWann*numWann; i++)
                fprintf(fev, "\t%.10lg\t%.10lg", std::real(rexp.wan[i]), std::imag(rexp.wan[i]));
            fprintf(fev, "\n");
        }
        if ( fclose(fev) ) {
            Logger::error("Expectation value file 'rho_%s' not written properly\n", rhoName.c_str());
            return false;
        }
    }
    return true;
}



bool saveSBEresults_txt(const ExpectationValues_t &ev, const Parameter_t &param, const std::string &subDirectory)
{
    namespace fs = std::filesystem;
    fs::path p = std::filesystem::path(param.out.saveDir);
    if ( subDirectory.size() > 0){
        auto np = p / fs::path(subDirectory);
        if ( ! fs::create_directory(np)){
            Logger::error("Could not create subdirectory '%s'\n", subDirectory.c_str());
            return false;
        }
        p = np;
    }

    unsigned numWann = param.tb.numWann;
    FILE * fPulse = fopen( ( p / "pulse.txt").c_str(), "w");
    if ( ! fPulse ){
        Logger::error("Could not open pulse.txt\n");
        return false;
    }
    fprintf(fPulse, "# all in atomic units\n");
    fprintf(fPulse, "# t\tEx\tEy\tEz\tAx\tAy\tAz\n");
    for(unsigned i=0; i<ev.t.size(); i++){
        fprintf(fPulse, "%.5lg", ev.t[i]);
        for(unsigned dir=0; dir<3; dir++)
            fprintf(fPulse, "\t%.10lg", ev.E[i][dir]);
        for(unsigned dir=0; dir<3; dir++)
            fprintf(fPulse, "\t%.10lg", ev.A[i][dir]);
        fprintf(fPulse, "\n");
    }
    if ( fclose(fPulse) ){
        Logger::error("Pulse file not written properly\n");
        return false;
    }

    FILE * fOccLimits = fopen( ( p / "occupationLimits.txt").c_str(), "w");
    if ( ! fPulse ){
        Logger::error("Could not open occupationLimits.txt\n");
        return false;
    }
    fprintf(fOccLimits, "# minimum and maximum occupation over k-points\n");
    fprintf(fOccLimits, "# t");
    for(unsigned bi=0; bi<numWann; bi++){
        fprintf(fOccLimits, "\toccHamMin%u", bi);
        fprintf(fOccLimits, "\toccHamMax%u", bi);
    }
    fprintf(fOccLimits, "\n");
    for(unsigned i=0; i<ev.t.size(); i++){
        fprintf(fOccLimits, "%lg", ev.t[i]);
        for(unsigned bi=0; bi<numWann; bi++){
            fprintf(fOccLimits, "\t%.10lg", ev.occupationHamMin[numWann*i+bi]);
            fprintf(fOccLimits, "\t%.10lg", ev.occupationHamMax[numWann*i+bi]);
        }
        fprintf(fOccLimits, "\n");
    }
    if ( fclose(fOccLimits) ){
        Logger::error("Occupation limit file not written properly\n");
        return false;
    }
    for(unsigned rid=0; rid<ev.region.size(); rid++){
        unsigned sampleIndex = 0;
        const std::array<unsigned, 3> dim = ev.region[rid].dim;
        for(unsigned a=0; a<dim[0]; a++)
            for(unsigned b=0; b<dim[1]; b++)
                for(unsigned c=0; c<dim[2]; c++){
                    std::string name = ev.region[rid].name;
                    if( dim[0] > 1)
                        name += "_" + std::to_string(a);
                    if( dim[1] > 1)
                        name += "_" + std::to_string(b);
                    if( dim[2] > 1)
                        name += "_" + std::to_string(c);
                    if ( ! saveRegionSample(param, p, name + ".txt", ev, rid, sampleIndex++) )
                        return false;
                }
    }
    return true;
}

std::vector<size_t> squeezeShape(const std::initializer_list<size_t> &l)
{
    std::vector<size_t> shape;
    for(const size_t v : l)
        if ( v > 1 )
            shape.push_back(v);
    return shape;
}


bool saveSBEresults_npz(const ExpectationValues_t &ev, const Parameter_t &param, std::string name)
{
    if (name.size() == 0){
        name = "data";
    }
    std::string fname = std::filesystem::path(param.out.saveDir) / (name + ".npz");
    size_t numWann = param.tb.numWann;
    size_t Nt = ev.t.size();
    size_t Ntdi = ev.tDensIndices.size();
    cnpy::npz_save(fname, "t", ev.t.data(), {Nt}, "w");
    cnpy::npz_save(fname, "pulse/E", ev.E[0].data(), {Nt, 3}, "a");
    cnpy::npz_save(fname, "pulse/A", ev.A[0].data(), {Nt, 3}, "a");
    cnpy::npz_save(fname, "lattice", unravel(param.tb.latticeVectors).data(), {3, 3}, "a");
    cnpy::npz_save(fname, "occHamMin", ev.occupationHamMin.data(), {Nt, numWann}, "a");
    cnpy::npz_save(fname, "occHamMax", ev.occupationHamMax.data(), {Nt, numWann}, "a");
    if ( param.ksr.storeDensityMatrix )
        cnpy::npz_save(fname, "tDensIndices", ev.tDensIndices.data(), {Ntdi}, "a");
    for(unsigned rid=0; rid<ev.region.size(); rid++){
        const RegionExpectationValues_t &rexp = ev.region[rid];
        const std::string &n = rexp.name;
        const size_t d[3] = {rexp.dim[0], rexp.dim[1], rexp.dim[2]};
        cnpy::npz_save(fname, n + "/weight", rexp.kPointWeightSum.data(), squeezeShape({Nt, d[0], d[1], d[2]}), "a");
        cnpy::npz_save(fname, n + "/j", rexp.j[0].data(), squeezeShape({Nt, d[0], d[1], d[2], 3}), "a");
        cnpy::npz_save(fname, n + "/occW", rexp.occupationWan.data(), squeezeShape({Nt, d[0], d[1], d[2], numWann}), "a");
        cnpy::npz_save(fname, n + "/occH", rexp.occupationHam.data(), squeezeShape({Nt, d[0], d[1], d[2], numWann}), "a");
        if ( param.ksr.storeDensityMatrix )
            cnpy::npz_save(fname, n + "/wan", rexp.wan.data(), squeezeShape({Ntdi, d[0], d[1], d[2], numWann, numWann}), "a");
    }
    return true;
}

bool saveSBEresults(const ExpectationValues_t &ev, const Parameter_t &param, const std::string &name)
{
    if ( param.out.saveAsNpz ) {
        return saveSBEresults_npz(ev, param, name);
    } else {
        return saveSBEresults_txt(ev, param, name);
    }
}
