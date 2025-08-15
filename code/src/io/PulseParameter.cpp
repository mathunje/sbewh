#include "io/PulseParameter.h"

void PulseParameter_t::print(FILE *f) const
{
    pulse.print(f);
}


/***********************************************/

PulseParameterValidator::PulseParameterValidator(ParameterIO &paramIO, const std::string &fixedInputDir) : fixedInputDir(fixedInputDir)
{
    std::function<bool(std::string)> callback = std::bind( std::mem_fn(&PulseParameterValidator::newMultiConfig), this, std::placeholders::_1);
    paramIO.startGroup("Pulse", "Defines a pulse, optionally as sum of multiple subpulses.\n"
                                "The default values for all subpulses may be set using the variables without number.\n"
                                "If no variable with pulse number is set, the default pulse will be used.\n"
                                "Otherwise each (partially) set subpulse will be used in the overall pulse.\n"
                                "The pulse run mode can be used to obtain the pulse without any SBE calculation.\n"
                                "The propagation time is defined by the earliest and latest pulse.\n"
                                "The maximum number of subpulses may be adapted in sbeConfigWh.h.in");

    configParam = new ParameterMultiConfig("multiRunFile",
                "This file must contain pulse specifications, which are seperated by [runName].\n"
                "\tThe global input file parameters are overridden by the specific runs.\n"
                "\tIf file name is not set or not existing, the global input file specifications will be used",
                 callback, fixedInputDir);
    paramIO.registerParam(configParam);


    paramIO.registerParam(new ParameterVector<double>("pol", "default polarization", pol.defaultValue, ", "));
    paramIO.registerParam(new ParameterString("type", "default pulse type: sin2, sin2Ramp, gauss, pureGauss", type.defaultValue));
    auto * pEmax = new ParameterNumeric<double>("Emax", "default maximum field strength in V/nm", Emax.defaultValue);
    pEmax->setMinValue(0);
    paramIO.registerParam(pEmax);
    auto * pLambda = new ParameterNumeric<double>("lambda", "default wavelength in nm", lambda.defaultValue);
    pLambda->setMinValue(0);
    paramIO.registerParam(pLambda);
    paramIO.registerParam(new ParameterNumeric<double>("cep", "default carrier envelope phase in radians", cep.defaultValue, 0));
    auto *pFWHM = new ParameterNumeric<double>("FWHM", "default full width half maximum of pulse in fs (except for sin^2 ramp pulse)",
                                                       tFWHM_fs.defaultValue);
    pFWHM->setMinValue(0);
    paramIO.registerParam(pFWHM);
    paramIO.registerParam(new ParameterNumeric<double>("tCentral", "default central time for single pulse in fs", tCentral_fs.defaultValue, 0));
    auto *pFWHMstart = new ParameterNumeric<double>("tStart", "Latest start time in units of default FWHM.\n"
                                                              "\tThe pulses are symmetrically defined around tCentral",
                                                              tStartFWHM.defaultValue);
    //pFWHMstart->setMaxValue(0);
    paramIO.registerParam(pFWHMstart);

    auto *pFWHMend = new ParameterNumeric<double>("tEnd", "Earliest end time in units of default FWHM.\n"
                                                          "\tThe pulses are symmetrically defined around tCentral",
                                                          tEndFWHM.defaultValue);
    // pFWHMend->setMinValue(0);
    paramIO.registerParam(pFWHMend);
    auto * pCo = new ParameterNumeric<int>("cyclesOn", "default number of cycles with maximum intensity in sin^2 ramp pulse",
                                                       cyclesOn.defaultValue, 3);
    pCo->setMinValue(1);
    paramIO.registerParam(pCo);
    auto * pRc = new ParameterNumeric<double>("risingCycles", "default number of rising and falling cycles in sin^2 ramp pulse",
                                                              risingCycles.defaultValue, 1);
    pRc->setMinValue(0);
    paramIO.registerParam(pRc);

    for(unsigned i=0; i<SBE_WH_MAX_PULSE_COUNT; i++){
        std::string nr = std::to_string(i);
        std::string des = "for pulse " + nr;
        auto * pPol = new ParameterVector<double>("pol"+nr, des, pol.pulseValue[i], ", ");
        pol.setPtr[i] = pPol->getSetPtr();
        paramIO.registerParam(pPol);

        auto * pType = new ParameterString("type"+nr, des, type.pulseValue[i], "");
        type.setPtr[i] = pType->getSetPtr();
        paramIO.registerParam(pType);

        auto * pEmax = new ParameterNumeric<double>("Emax"+nr, des, Emax.pulseValue[i], 0);
        pEmax->setMinValue(0);
        Emax.setPtr[i] = pEmax->getSetPtr();
        paramIO.registerParam(pEmax);

        auto * pLambda = new ParameterNumeric<double>("lambda"+nr, des, lambda.pulseValue[i], 0);
        pLambda->setMinValue(0);
        lambda.setPtr[i] = pLambda->getSetPtr();
        paramIO.registerParam(pLambda);

        auto *pCep = new ParameterNumeric<double>("cep"+nr, des, cep.pulseValue[i], 0);
        cep.setPtr[i] = pCep->getSetPtr();
        paramIO.registerParam(pCep);

        auto *pFWHM = new ParameterNumeric<double>("FWHM"+nr, des, tFWHM_fs.pulseValue[i], 0);
        //pFWHM->setMinValue(0);
        tFWHM_fs.setPtr[i] = pFWHM->getSetPtr();
        paramIO.registerParam(pFWHM);

        auto * pCentral = new ParameterNumeric<double>("tCentral"+nr, des, tCentral_fs.pulseValue[i], 0);
        tCentral_fs.setPtr[i] = pCentral->getSetPtr();
        paramIO.registerParam(pCentral);

        auto *pFWHMstart = new ParameterNumeric<double>("tStart"+nr, des, tStartFWHM.pulseValue[i], 0);
        //pFWHMstart->setMaxValue(0);
        tStartFWHM.setPtr[i] = pFWHMstart->getSetPtr();
        paramIO.registerParam(pFWHMstart);

        auto *pFWHMend = new ParameterNumeric<double>("tEnd"+nr, des, tEndFWHM.pulseValue[i], 0);
        // pFWHMend->setMinValue(0);
        tEndFWHM.setPtr[i] = pFWHMend->getSetPtr();
        paramIO.registerParam(pFWHMend);

        auto * pCo = new ParameterNumeric<int>("cyclesOn"+nr, des, cyclesOn.pulseValue[i], 3);
        pCo->setMinValue(1);
        cyclesOn.setPtr[i] = pCo->getSetPtr();
        paramIO.registerParam(pCo);

        auto * pRc = new ParameterNumeric<double>("risingCycles"+nr, des, risingCycles.pulseValue[i], 1);
        pRc->setMinValue(0);
        risingCycles.setPtr[i] = pRc->getSetPtr();
        paramIO.registerParam(pRc);
    }
    paramIO.finishGroup();
}

bool PulseParameterValidator::addSubPulse(PulseParameter_t *p, unsigned dim, unsigned nr)
{
    const auto &vPol = pol.get(nr);
    if ( vPol.size() != 3 ){
        Logger::error("Polarization for pulse %d\n must have exactly 3 entries (not %lu)\n", (int)nr, vPol.size() );
        return false;
    }
    GeomVector3d cpol({vPol[0], vPol[1], vPol[2]});
    for(unsigned i=dim; i<3; i++)
        if ( std::abs(vPol[i]) > 1e-10 ){
            Logger::warn("Polarization of subpulse %d has non-zero component in %uth direction\n", (int)nr, dim);
        }
    std::string ctype = type.get(nr);
    double cEmax = atomicUnits::from_V_nm( Emax.get(nr) );
    double omega = 2 * M_PI * atomicUnits::speedOfLight / atomicUnits::from_nm( lambda.get(nr) );
    double ccep = cep.get(nr);
    double FWHM = atomicUnits::from_fs( tFWHM_fs.get(nr) );
    double central = atomicUnits::from_fs ( tCentral_fs.get(nr) );
    double startFWHM = tStartFWHM.get(nr);
    double endFWHM = tEndFWHM.get(nr);
    if ( startFWHM >= endFWHM ){
        Logger::error("tStart must be before tEnd (pulse %d)\n", (int)nr);
        return false;
    }
    int ccyclesOn = cyclesOn.get(nr);
    double crisingCycles = risingCycles.get(nr);
    std::shared_ptr<Pulse1D> pulse1D = nullptr;
    if ( ! ctype.compare("gauss") )
        pulse1D = std::shared_ptr<Pulse1D>(new GaussianPulse1D(cEmax, omega, ccep, FWHM, startFWHM, endFWHM));
    if ( ! ctype.compare("pureGauss") )
        pulse1D = std::shared_ptr<Pulse1D>(new PureGauss1D(cEmax, FWHM, startFWHM, endFWHM));
    if ( ! ctype.compare("sin2") )
        pulse1D = std::shared_ptr<Pulse1D>(new Sin2Pulse1D(cEmax, omega, ccep, FWHM, startFWHM, endFWHM));
    if ( ! ctype.compare("sin2Ramp") )
        pulse1D = std::shared_ptr<Pulse1D>(new Sin2RampPulse1D(cEmax, omega, crisingCycles, ccyclesOn));
    if ( ! pulse1D ){
        Logger::error("Unrecognized pulse envelope type '%s' for pulse %d\n", ctype.c_str(), (int)nr);
        return false;
    }
    p->pulse.addPulse(cpol, pulse1D, central);
    return true;
}

bool PulseParameterValidator::update(PulseParameter_t *p, unsigned dim)
{
    p->multiRunFname = configParam->getFileName();
    p->pulse.clear();

    if ( pol.defaultValue.size() > 0 && pol.defaultValue.size() != 3 ){
        Logger::error("Default polarization must have 3 entries (not %lu)\n", pol.defaultValue.size());
        return false;
    }
    if ( tStartFWHM.defaultValue >= tEndFWHM.defaultValue ){
        Logger::error("Simulation start time has to be before simulation end time)");
        return false;
    }

    bool anyPulseSet = false;
    for(unsigned i=0; i<SBE_WH_MAX_PULSE_COUNT; i++){
        bool pulseSet = pol.isSet(i) || type.isSet(i)|| Emax.isSet(i) || lambda.isSet(i) || cep.isSet(i) ||
                        tFWHM_fs.isSet(i) || tCentral_fs.isSet(i) || tStartFWHM.isSet(i) || tEndFWHM.isSet(i) ||
                        cyclesOn.isSet(i) || risingCycles.isSet(i);
        if ( ! pulseSet )
            continue;
        anyPulseSet = true;
        if ( ! addSubPulse(p, dim, i) ){
            Logger::error("Definition for pulse %d is ill-formed\n", i);
            return false;
        }
    }
    if ( !anyPulseSet ){
        Logger::verbose("No sub pulse defined -- using default pulse parameters to construct pulse\n");
        if ( ! addSubPulse(p, dim, -1) ){
            Logger::error("Definition of default pulse is ill-formed\n");
            return false;
        }
    }
    return true;
}

bool PulseParameterValidator::newMultiConfig(std::string name)
{
    PulseParameter_t lp;
    bool valid = update(&lp);
    if ( valid ) {
        configurations.insert({name, lp});
    } else {
        configurations.clear();
    }
    return valid;
}

void mpi_bcast_Pulse1D(std::shared_ptr<Pulse1D> &lsr, const MpiParameter_t &mp){
    Pulse1Dtype type;
    if ( mp.rank == mp.root )
        type = lsr->getType();
    mpi_bcast_trivial(type, mp);
    switch(type){
        case Gaussian: {
                GaussianPulse1Draw_t raw;
                if ( mp.rank == mp.root )
                    raw = static_cast<GaussianPulse1D*>(lsr.get())->getRaw();
                mpi_bcast_trivial(raw, mp);
                if ( mp.rank != mp.root )
                    lsr = std::shared_ptr<Pulse1D>(new GaussianPulse1D(raw));
            }
            break;
        case PureGauss: {
                PureGauss1Draw_t raw;
                if ( mp.rank == mp.root )
                    raw = static_cast<PureGauss1D*>(lsr.get())->getRaw();
                mpi_bcast_trivial(raw, mp);
                if ( mp.rank != mp.root )
                    lsr = std::shared_ptr<Pulse1D>(new PureGauss1D(raw));
            }
            break;
        case Sin2Pulse:{
                Sin2Pulse1Draw_t raw;
                if ( mp.rank == mp.root )
                    raw = static_cast<Sin2Pulse1D*>(lsr.get())->getRaw();
                mpi_bcast_trivial(raw, mp);
                if ( mp.rank != mp.root )
                    lsr = std::shared_ptr<Pulse1D>(new Sin2Pulse1D(raw));
            }
            break;
        case Sin2RampPulse: {
                Sin2RampPulse1Draw_t raw;
                if ( mp.rank == mp.root )
                    raw = static_cast<Sin2RampPulse1D*>(lsr.get())->getRaw();
                mpi_bcast_trivial(raw, mp);
                if ( mp.rank != mp.root )
                    lsr = std::shared_ptr<Pulse1D>(new Sin2RampPulse1D(raw));
            }
            break;
    }
}

void mpi_bcast_PulseParameter(PulseParameter_t &p, const MpiParameter_t &mp)
{
    mpi_bcast(p.multiRunFname, mp);
    unsigned pulseCount = 0;
    if ( mp.rank == mp.root )
        pulseCount = p.pulse.getSinglePulseCount();
    mpi_bcast_trivial(pulseCount, mp);
    for(unsigned i=0; i<pulseCount; i++){
        SinglePulse3D_t sp;
        if ( mp.rank == mp.root )
            sp = p.pulse.getSinglePulse(i);
        mpi_bcast_trivial(sp.center, mp);
        mpi_bcast(sp.pol, mp);
        mpi_bcast_Pulse1D(sp.pulse1D, mp);
        if ( mp.rank != mp.root )
            p.pulse.addPulse(sp);
    }
}
