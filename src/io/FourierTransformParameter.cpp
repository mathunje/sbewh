#include "io/FourierTransformParameter.h"



void FourierTransformParameter_t::print(FILE *f, unsigned dim) const
{
    fprintf(f, "Fouriertransform parameter\n");
    fprintf(f, "\tperform direct calculation: %d\n", int(performDirectCalculation));
    fprintf(f, "\tHamiltonian with test field: %d\n", int(whTransformTestField));
    for(unsigned d=0; d<dim; d++)
        fprintf(f, "\tNf%u: %u\n", d+1, Nf[d]);
    std::string pm = FourierTransformParameterValidator::planningFlagToStr(fftwPlanningFlag);
    fprintf(f, "\tFFTW planing mode: %s\n", pm.c_str());
    fprintf(f, "\tFFTW planing time limit: %lg s\n", planningTimeLimit);
    fprintf(f, "\tsaveIndices:");
    for(auto s : saveIndices)
        fprintf(f, " %u", s);
    fprintf(f, "\n");
}


std::map<std::string, unsigned> FourierTransformParameterValidator::planningModeMap;

FourierTransformParameterValidator::FourierTransformParameterValidator(ParameterIO &paramIO)
{
    if ( ! planningModeMap.size() ){
        planningModeMap.insert({"estimate", FFTW_ESTIMATE});
        planningModeMap.insert({"measure", FFTW_MEASURE});
        planningModeMap.insert({"patient", FFTW_PATIENT});
        planningModeMap.insert({"exhaustive", FFTW_EXHAUSTIVE});
    }

    paramIO.startGroup("FourierTransform", "Control usage of FFTs for calculation of Hamiltonian and dipole matrix elements");
    paramIO.registerParam(new ParameterBool("performDirectCalculation",
                "if set, the Fourier transform is evalualuated via direct sum for comparison, too.", performDirectCalculation, false));
    paramIO.registerParam(new ParameterBool("whTransformTestField",
                "if set, the dipole matrix elements are obtained from full Hamiltonian via a test field\n"
                "\tin WHtransform mode. Mainly used for debugging.", whTransformTestField, false));
    paramIO.registerParam(new ParameterVector<unsigned>("saveIndices",
                "defines submatrix for which H and D are saved in WHtransform run mode (0-indexed).\n"
                "\tIf this parameter is not set, all matrix elements will be written\n"
                "\tWARNING: This may lead to huge files",  saveIndices, {}, ", "));
    std::string description = "planning mode used by fftw, must be one of:";
    for(const auto &[s, _] : planningModeMap)
        description += " " + s;
    paramIO.registerParam(new ParameterString("fftPlanningMode", description, fftPlanningMode, "measure"));
    paramIO.registerParam(new ParameterNumeric<double>("planningTimeLimit", "time limit in seconds for fftw to search for best plan",
                planningTimeLimit, 1));
    paramIO.registerParam(new ParameterVector<unsigned>("Nf",
                "defines number of k points in each direction (coarse submesh), for which FFT may be performed\n"
                "\tIf not set (or smaller than extent of Wannier file it is automatically set.", Nf, {1}, ", "));
    paramIO.registerParam(new ParameterNumeric<unsigned>("maxNfPrime",
                "each dimension of Nf is increased until all prime factors are at most the provided value", maxNfPrime, 5));
    paramIO.finishGroup();
}

bool FourierTransformParameterValidator::update(FourierTransformParameter_t *p, const TightBindingParameter_t &tb, unsigned dim)
{
    p->performDirectCalculation = performDirectCalculation;
    p->whTransformTestField = whTransformTestField;
    if ( saveIndices.size() == 0){
        for(unsigned i=0; i<tb.numWann; i++)
            saveIndices.push_back(i);
    }
    for(unsigned si : saveIndices)
        if ( si >= tb.numWann){
            Logger::error("To high saveIndex for FourierTransform %u/%u\n", si, tb.numWann);
            return false;
        }
    p->saveIndices = saveIndices;
    auto it = planningModeMap.find(fftPlanningMode);
    if ( it == planningModeMap.end()){
        Logger::error("Unrecognized fftw planning mode '%s'\n", fftPlanningMode.c_str());
        return false;
    }
    p->fftwPlanningFlag = it->second;
    p->planningTimeLimit = planningTimeLimit;

    p->maxNfPrime = maxNfPrime;
    if ( Nf.size() != 1 && Nf.size() < dim ){
        Logger::error("Wrong number (%llu) of dimensions in Nf\n", Nf.size());
        return false;
    }
    if ( Nf.size() == 1){
        for(unsigned i=0; i<dim; i++)
           p->Nf[i] = Nf[0];
    } else {
        for(unsigned i=0; i<dim; i++)
           p->Nf[i] = Nf[i];
        for(unsigned i=dim; i<Nf.size(); i++)
            if ( Nf[i] != 1){
                Logger::error("%d. component of Nf must be 1 or omitted (dimension: %d)\n", i+1, dim);
                return false;
            }
    }

    for(unsigned i=dim; i<3; i++)
        p->Nf[i] = 1;
    return true;
}


std::string FourierTransformParameterValidator::planningFlagToStr(unsigned flag)
{
    for(const auto &[s, pFlag] : planningModeMap)
        if ( flag == pFlag )
            return s;
    return "unkown";
}


void mpi_bcast_FourierTransformParameter(FourierTransformParameter_t &p, const MpiParameter_t &mp)
{
    mpi_bcast_trivial(p.performDirectCalculation, mp);
    mpi_bcast_trivial(p.whTransformTestField, mp);
    mpi_bcast_trivial(p.saveIndices, mp);
    mpi_bcast_trivial(p.Nf, mp);
    mpi_bcast_trivial(p.maxNfPrime, mp);
    mpi_bcast_trivial(p.fftwPlanningFlag, mp);
    mpi_bcast_trivial(p.planningTimeLimit, mp);
}
