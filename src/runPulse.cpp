#include "runModes.h"
#include "Parameter.h"
#include "cnpy.h"

#include<vector>

#include<boost/numeric/odeint.hpp>
#include "GeomVector.hpp"


RequiredParameter_t getRequiredParameterPulse()
{
    return RequiredParameter_t {
              .par = true,
              .tb = false,
              .ft = false,
              .diag = false,
              .pulse = true,
              .prop = true,
              .ksr = false,
              .out = true,
              .multiConfigs = false
           };
}

typedef std::vector<double> PulseState;

struct PulseObserver {
    unsigned dim;
    const Pulse3D & pulse;
    std::vector<double> E;
    std::vector<double> A;
    std::vector<double> time;

    PulseObserver(unsigned dim, const Pulse3D & pulse) : dim(dim), pulse(pulse) {}
    void operator()(const PulseState &s, double t)
    {
        time.push_back(t);
        for(unsigned i=0; i<pulse.getSinglePulseCount(); i++){
            const GeomVector3d Et = pulse.getSinglePulseE(i, t);
            for(unsigned dir=0; dir<dim; dir++)
                E.push_back(Et.at(dir));
        }
        for(double At : s)
            A.push_back(At);

    }
};

struct PulseDeriv {
    unsigned dim;
    const Pulse3D &pulse;
    PulseDeriv(unsigned dim, const Pulse3D &pulse) : dim(dim), pulse(pulse) {}

    PulseState getInitialState() { return PulseState(dim * pulse.getSinglePulseCount(), 0); }

    void operator() (const PulseState &state, PulseState &dsdt, const double t)
    {
        (void)(state);
        for(unsigned i=0; i<pulse.getSinglePulseCount(); i++){
            const GeomVector3d E = pulse.getSinglePulseE(i, t);
            for(unsigned dir=0; dir<dim; dir++)
                dsdt[dim*i+dir] = -E.at(dir);
        }
    }
};

int runPulse(Parameter_t &param){
    unsigned dim = param.general.dim;
    Pulse3D combinedPulse = param.pulse.pulse;
    PulseDeriv pd(dim, combinedPulse);
    PulseObserver observer(dim, combinedPulse);

    using namespace boost::numeric::odeint;
    PulseState initialState = pd.getInitialState();
    const PropagationParameter_t & ppp = param.prop;
    double timeStep = (ppp.endTime - ppp.startTime) / (ppp.outputCount - 1);
    integrate_const(
                    make_controlled < runge_kutta_dopri5< PulseState > >(ppp.epsAbs, ppp.epsRel),
                    std::ref(pd), initialState, ppp.startTime, ppp.endTime + timeStep/2, timeStep, std::ref(observer) );

    unsigned size = observer.time.size();
    unsigned pulseCount =  combinedPulse.getSinglePulseCount();
    std::filesystem::path p(param.out.saveDir);
    if ( param.out.saveAsNpz ){
        std::string fname = std::string(p / "pulseDetail.npz");
        cnpy::npz_save(fname, "t", &observer.time[0], {size}, "w");
        cnpy::npz_save(fname, "E", &observer.E[0], {size, pulseCount, dim}, "a");
        cnpy::npz_save(fname, "A", &observer.A[0], {size, pulseCount, dim}, "a");
    } else {
        std::string fname = std::string(p / "pulseDetail.txt");
        FILE * f = fopen(fname.c_str(), "w");
        if ( ! f ){
            Logger::error("Could not open '%s'\n", fname.c_str());
            return 1;
        }
        fprintf(f, "# Subpulses in atomic units\n");
        fprintf(f, "# t");
        const char dirs[3] = {'x', 'y', 'z'};
        for(unsigned i=0; i<pulseCount; i++)
            for(unsigned d=0; d<dim; d++)
                fprintf(f, "\tE%u%c", i, dirs[d]);
        for(unsigned i=0; i<pulseCount; i++)
            for(unsigned d=0; d<dim; d++)
                fprintf(f, "\tA%u%c", i, dirs[d]);
        fprintf(f, "\n");
        for(unsigned t=0; t<size; t++){
            fprintf(f, "%.5lg\n", observer.time[t]);
            for(unsigned i=0; i<pulseCount; i++)
                for(unsigned d=0; d<dim; d++)
                    fprintf(f, "\t%.10lg", observer.E[t*pulseCount*dim+i*dim+d]);
            for(unsigned i=0; i<pulseCount; i++)
                for(unsigned d=0; d<dim; d++)
                    fprintf(f, "\t%.10lg", observer.A[t*pulseCount*dim+i*dim+d]);
            fprintf(f, "\n");
        }
        if ( fclose(f) ){
            Logger::error("Could not write '%s'\n", fname.c_str());
            return 1;
        }
    }
    return 0;
}
