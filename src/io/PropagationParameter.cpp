#include "io/PropagationParameter.h"

void PropagationParameter_t::print(FILE *f)
{
    fprintf(f, "Propagation parameter\n");
    fprintf(f, "\toutputCount: %u\n", outputCount);
    fprintf(f, "\tstartTime: %lf fs\n", atomicUnits::to_fs(startTime));
    fprintf(f, "\tendTime: %lf fs\n", atomicUnits::to_fs(endTime));
    fprintf(f, "\tebsAbs: %g\n", epsAbs);
    fprintf(f, "\tebsRel: %g\n", epsRel);
    for(unsigned i=0; i<3; i++)
        fprintf(f, "\tNk%u: %d\n", i, Nk[i]);
    fprintf(f, "\tT2: %lg fs\n", atomicUnits::to_fs(1/relaxationF2));
    fprintf(f, "\tsoothingWidth: %lg eV\n", atomicUnits::to_eV(soothingWidth));
    fprintf(f, "\toccupationSmearingWidth: %lg eV\n", atomicUnits::to_eV(occupationSmearingWidth));
    fprintf(f, "\tkGlobalShift: %.4lg %.4lg %.4lg\n", kGlobalShift[0], kGlobalShift[1], kGlobalShift[2]);
    fprintf(f, "\n");
    double pi = std::numbers::pi;
    double omegaMax = pi * (outputCount-1) / (endTime - startTime);
    fprintf(f, "\tDerived Nyquist frequency: %.4lg (%.2lg eV)\n", omegaMax, atomicUnits::to_eV(omegaMax));
    double omegaRes = 2 * pi / (endTime - startTime);
    fprintf(f, "\tEstimated resolution: %.6lg (%.9lg eV)\n", omegaRes, atomicUnits::to_eV(omegaRes));
}

/*****************************************/


PropagationParameterValidator::PropagationParameterValidator(ParameterIO &paramIO)
{
    paramIO.startGroup("Propagation", "for numerical integration of SBEs\n");
    paramIO.registerParam(new ParameterNumeric<unsigned>("outputCount", "number of time steps in output", outputCount, 500));
    paramIO.registerParam(new ParameterVector<unsigned>("Nk", "array of number of simulation grid points in each k direction", Nk, {1, 1, 1}, ", "));
    paramIO.registerParam(new ParameterBool("allowNkAdaption", "if allowed: adjusts the Nks to be multiple of Nf\n"
                                                               "\totherwise error is reported", allowNkAdaption, false));
    paramIO.registerParam(new ParameterVector<double>("kGlobalShift", "global shift of k-space Fourier meshes in fractional coordinates", kGlobalShift,
                          {0, 0, 0}, ", "));

    auto pea = new ParameterNumeric<double>("epsAbs",
            "absolute error bound for numerical integration per step", epsAbs, 1e-10);
    pea->setRange(0, 1e-2);
    paramIO.registerParam(pea);
    auto per = new ParameterNumeric<double>("epsRel",
            "relative error bound for numerical integration per step", epsRel, 1e-10);
    per->setRange(0, 1e-2);
    paramIO.registerParam(per);

    auto pt2 = new ParameterNumeric<double>("T2", "dephasing time constant (fs) - a negative value means no dephasing at all", T2_fs, -1);
    paramIO.registerParam(pt2);
    auto psw = new ParameterNumeric<double>("soothingWidth", "gaussian witdh of reduction of dephasing (in energy differences) of the dephasing operator\n"
                                            "\ta negative value indicates, that soothing is disabled",
                                            soothingWidth_eV, 0.025);
    psw->setRange(-10, 10);
    paramIO.registerParam(psw);
    auto posw = new ParameterNumeric<double>("occupationSmearingWidth", "width of gaussian like weighting applied to Hamiltonian occupation numbers (eV)\n"
                                            "\tdefaults to soothingWidth\n"
                                            "\ta negative value indicates, that smearing is disabled",
                                            occupationSmearingWidth_eV, 0.025);
    posw->setRange(-10, 10);
    occupationSmearingWidthSetPtr = posw->getSetPtr();
    paramIO.registerParam(posw);
    paramIO.finishGroup();
}


bool PropagationParameterValidator::update(PropagationParameter_t *p, const PulseParameter_t &pulse, unsigned dim)
{
    p->outputCount = outputCount;
    if ( Nk.size() != 1 &&  Nk.size() < dim){
        Logger::error("Wrong number (%llu) of dimensions in Nk\n", Nk.size());
        return false;
    }
    if ( kGlobalShift.size() > 3){
        Logger::error("kGlobalShift must have at most three dimensions\n");
        return false;
    }
    for(unsigned dir=0; dir<kGlobalShift.size(); dir++){
        double ck = kGlobalShift[dir];
        if ( ck < -1e-10 || ck > 1 + 1e-10){
            Logger::error("kGlobalShift entry %d (= %.4lg) must be in range [0, 1]\n", dir+1, ck);
            return false;
        }
        p->kGlobalShift[dir] = ck;
    }
    if ( Nk.size() == 1){
        for(unsigned i=0; i<dim; i++)
            p->Nk[i] = Nk[0];
    } else {
         for(unsigned i=0; i<dim; i++)
            p->Nk[i] = Nk[i];
        for(unsigned i=dim; i<Nk.size(); i++)
            if ( Nk[i] != 1){
                Logger::error("%d. component of Nk must be 1 or omitted (dimension: %d)\n", i+1, dim);
                return false;
            }
    }
    for(unsigned i=dim; i<3; i++)
        p->Nk[i] = 1;
    p->allowNkAdaption = allowNkAdaption;
    p->startTime = pulse.pulse.getStartTime();
    p->endTime = pulse.pulse.getEndTime();
    p->epsAbs = epsAbs;
    p->epsRel = epsRel;
    p->relaxationF2 = 1 / atomicUnits::from_fs(T2_fs);
    p->soothingWidth = atomicUnits::from_eV(soothingWidth_eV);
    if ( *occupationSmearingWidthSetPtr ){
        p->occupationSmearingWidth = atomicUnits::from_eV(occupationSmearingWidth_eV);
    } else {
        p->occupationSmearingWidth = p->soothingWidth;
    }
    return true;
}

void mpi_bcast_PropagationParameter(PropagationParameter_t &p, const MpiParameter_t &mp)
{
    mpi_bcast_trivial(p.outputCount, mp);
    mpi_bcast_trivial(p.startTime, mp);
    mpi_bcast_trivial(p.endTime, mp);
    mpi_bcast_trivial(p.epsAbs, mp);
    mpi_bcast_trivial(p.epsRel, mp);
    mpi_bcast_trivial(p.relaxationF2, mp);
    mpi_bcast_trivial(p.soothingWidth, mp);
    mpi_bcast_trivial(p.occupationSmearingWidth, mp);
    mpi_bcast_trivial(p.Nk, mp);
    mpi_bcast_trivial(p.Noffset, mp);
    mpi_bcast(p.kGlobalShift, mp);
    mpi_bcast_trivial(p.allowNkAdaption, mp);
}
