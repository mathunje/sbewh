#include "Parameter.h"


bool validateKdimensions(Parameter_t &param)
{
    for(unsigned d=0; d<param.general.dim; d++){
        if ( param.ft.Nf[d] < param.tb.Nw[d] )
            param.ft.Nf[d] = param.tb.Nw[d];
    }
    // increase until Nf contains only small prime divisors
    unsigned maxNf = *std::max_element(param.ft.Nf.begin(), param.ft.Nf.end());
    if ( maxNf > param.ft.maxNfPrime){
        std::vector<unsigned> primes({2});
        for(unsigned s=3; s<=param.ft.maxNfPrime; s++){
            unsigned d = 2;
            while (d < s && (s%d != 0))
                d++;
            if ( d == s )
                primes.push_back(s);
        }
        for(unsigned d=0; d<param.general.dim; d++){
            unsigned N = param.ft.Nf[d] - 1;
            unsigned Nc = 1;
            do {
                Nc = ++N;
                for(auto p : primes)
                    while ( Nc % p == 0)
                        Nc /= p;
            } while( Nc != 1);
            param.ft.Nf[d] = N;
        }
    }
    for(unsigned d=0; d<param.general.dim; d++){
        if ( param.ft.Nf[d] == 1)
            continue;
        if ( param.prop.Nk[d] < param.ft.Nf[d]){
            Logger::error("In k-dimension %u: Nk (%u) is smaller than Nf (%u)\n",
                          d+1, param.prop.Nk[d], param.ft.Nf[d]);
            return false;
        }
        unsigned remainder = param.prop.Nk[d] % param.ft.Nf[d];
        if ( remainder != 0){
            if ( ! param.prop.allowNkAdaption ){
                Logger::error("In k-dimension %u: Nk (%u) is not a multiple of Nf (%u)\n",
                              d+1, param.prop.Nk[d], param.ft.Nf[d]);
                return false;
            }
            param.prop.Nk[d] -= remainder;
            if ( remainder >= param.ft.Nf[d]/2 )
                param.prop.Nk[d] += param.ft.Nf[d];
        }
    }
    for(unsigned d=0; d<3; d++)
        param.prop.Noffset[d] = param.prop.Nk[d] / param.ft.Nf[d];
    return true;
}

ParameterIO::ParseResult parseParameter(int argc, const char * argv[], Parameter_t &param, const std::string fixedInputDir)
{
    ParameterIO paramIO;
    GeneralParameterValidator genValidator(paramIO);
    ParallelParameterValidator parValidator(paramIO);
    PulseParameterValidator pulseValidator(paramIO, fixedInputDir);
    FourierTransformParameterValidator ftValidator(paramIO);
    DiagonalizationParameterValidator diagValidator(paramIO);
    TightBindingParameterValidator tbValidator(paramIO, fixedInputDir);
    PropagationParameterValidator propValidator(paramIO);
    KspaceRegionParameterValidator ksrValidator(paramIO);
    OutputParameterValidator outputValidator(paramIO);

    constexpr ParameterIO::ParseResult parseOk = ParameterIO::ParseResult::Ok;
    constexpr ParameterIO::ParseResult parseFail = ParameterIO::ParseResult::Fail;
    auto pr = paramIO.parse(argc, argv, "input.txt");
    if ( pr != parseOk ){
        if ( pr == parseFail )
            paramIO.printDiagnose(stdout);
        return pr;
    }
    param.cmdInp.fnames = paramIO.getParamFileNames();
    for(int i=1; i<argc; i++)
        param.cmdInp.args.push_back(argv[i]);

    if ( ! genValidator.update(&param.general) ){
        Logger::error("Defined general options are inconsistent\n");
        return parseFail;
    }
    if ( param.general.verbose )
        Logger::setVerbosity(Logger::Severity::VERBOSE);
    if ( param.general.veryVerbose )
        Logger::setVerbosity(Logger::Severity::VERYVERBOSE);

    RequiredParameter_t parse = getRequiredParameter(param.general.runMode);
    if( parse.out && ! outputValidator.update(&param.out, param.general.finishMultiRuns) ){
        Logger::error("Defined output options are inconsistent\n");
        return parseFail;
    }

    if ( param.general.finishMultiRuns ){
        OutputParameter_t outOrig = param.out;
        bool ok;
        std::vector<std::string> argsData = loadArgsFromInputDirectory(param.out.saveDir, &ok);
        if ( ! ok )
            return parseFail;
        argsData.insert(argsData.begin(), argv[0]);
        std::vector<const char*> args(argsData.size());
        const char ** lastArgv = args.data();
        for(unsigned i=0; i<args.size(); i++)
            lastArgv[i] = argsData[i].c_str();
        param.cmdInp.args.clear();
        pr = parseParameter((int)args.size(), lastArgv, param, param.out.saveDir);
        if ( pr != parseOk)
            return pr;
        if ( param.general.finishMultiRuns ){
            Logger::error("Recursive call of finishMultiRuns is not allowed\n");
            return parseFail;
        }
        param.out = outOrig;
        param.general.finishMultiRuns = true;
    }

    if ( parse.par && ! parValidator.update(&param.par, param.mpi) ){
        Logger::error("Defined parallelization options can not be fullfilled\n");
        return parseFail;
    }

    if ( parse.tb && ! tbValidator.update(&param.tb, param.general.dim) ){
        Logger::error("Defined tight binding options are inconsistent\n");
        return parseFail;
    }

    if ( parse.ft && ! ftValidator.update(&param.ft, param.tb, param.general.dim) ){
        Logger::error("Defined Fourier transform options are inconsistent\n");
        return parseFail;
    }

    if ( parse.diag && ! diagValidator.update(&param.diag) ) {
        Logger::error("Defined diagonalization options are inconsistent\n");
        return parseFail;
    }

    if ( parse.pulse && ! pulseValidator.update(&param.pulse, param.general.dim) ){
        Logger::error("Defined pulse paramters are inconsistent\n");
        return parseFail;
    }

    if ( parse.prop && ! propValidator.update(&param.prop, param.pulse, param.general.dim) ){
        Logger::error("Defined propagation parameters are inconsistent\n");
        return parseFail;
    }

    if ( parse.ksr && ! ksrValidator.update(&param.ksr) ){
        Logger::error("Defined k-space region parameters are inconsistent\n");
        return parseFail;
    }

    if ( parse.multiConfigs && ! paramIO.parseMultiConfigs() ){
        Logger::error("A configuration is invalid\n");
        return parseFail;
    }

    if ( parse.multiConfigs ){
        param.pulseConfigurations = pulseValidator.getConfigurations();
        for(const auto &[name, l] : param.pulseConfigurations){
            if ( param.general.veryVerbose ){
                Logger::print("\nConfiguration: %s\n", name.c_str());
                l.print(stdout);
            }
        }
    }

    if ( parse.tb && parse.ft && parse.prop && ! validateKdimensions(param) ){
        Logger::error("Defined dimensions in k-Space are not supported\n");
        return parseFail;
    }

    if ( param.general.veryVerbose ){
        printf("------------------------------ raw parsed input ------------------------------\n");
        paramIO.printParam(stdout);
        printf("------------------------------------------------------------------------------\n");
    }
    if ( param.general.verbose){
        printf("------------------------------- input params ---------------------------------\n");
        if ( parse.ft )
            param.ft.print(stdout, param.general.dim);
        if ( parse.tb )
            param.tb.print(stdout, param.general.dim);
        if ( parse.prop )
            param.prop.print(stdout);
        if ( parse.pulse )
            param.pulse.print(stdout);
        if ( parse.ksr )
            param.ksr.print(stdout);
        printf("------------------------------------------------------------------------------\n");
    }
    return parseOk;
}

void mpi_bcast_Params(Parameter_t &param)
{
    mpi_bcast(param.cmdInp.fnames, param.mpi);
    mpi_bcast(param.cmdInp.args, param.mpi);
    mpi_bcast_GeneralParameter(param.general, param.mpi);
    RequiredParameter_t bcast = getRequiredParameter(param.general.runMode);
    if ( bcast.par )
        mpi_bcast_ParallelParameter(param.par, param.mpi);
    if ( bcast.tb )
        mpi_bcast_TightBindingParameter(param.tb, param.mpi);
    if ( bcast.ft )
        mpi_bcast_FourierTransformParameter(param.ft, param.mpi);
    if ( bcast.diag )
        mpi_bcast_DiagonalizationParameter(param.diag, param.mpi);
    if ( bcast.prop )
        mpi_bcast_PropagationParameter(param.prop, param.mpi);
    if ( bcast.pulse )
        mpi_bcast_PulseParameter(param.pulse, param.mpi);
    if ( bcast.ksr )
        mpi_bcast_KspaceRegionParameter(param.ksr, param.mpi);
    if ( bcast.multiConfigs ){
        unsigned size = (unsigned)param.pulseConfigurations.size();
        mpi_bcast_trivial(size, param.mpi);
        auto it = param.pulseConfigurations.begin();
        for(unsigned i=0; i<size; i++){
            PulseParameter_t lp;
            std::string name;
            if(param.mpi.rank == param.mpi.root ){
                name = it->first;
                lp = it->second;
                it++;
            }
            mpi_bcast_PulseParameter(lp, param.mpi);
            mpi_bcast(name, param.mpi);
            if ( param.mpi.rank != param.mpi.root )
               param.pulseConfigurations.insert({name, lp});
        }
    }
    if ( bcast.out )
        mpi_bcast_OutputParameter(param.out, param.mpi);
}
