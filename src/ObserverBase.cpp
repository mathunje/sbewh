#include "ObserverBase.h"

ObserverBase::ObserverBase(const StateContext &sc) : sc(sc)
{
    for(const KspaceRegion_t &r : sc.param.ksr.regions){
        for(unsigned a=0; a<r.sampleDim[0]; a++)
            for(unsigned b=0; b<r.sampleDim[1]; b++)
                for(unsigned c=0; c<r.sampleDim[2]; c++){
                    GeomVector3d k = r.kFrac + (double)(a) * r.dkFrac[0] + (double)(b) * r.dkFrac[1] + (double)(c) * r.dkFrac[2];
                    region.push_back( {.k = { k[0], k[1], k[2]},
                                       .movingFrame = r.movingFrame,
                                       .invSigma2 = 8 * std::log(2) / (r.fwhm * r.fwhm)  });
                }
        regionContext.push_back({.name = r.name,
                                 .dim = r.sampleDim});
    }
}


void ObserverBase::basicLog(double t)
{
    if( sc.param.general.verbose && ! sc.param.par.independentMpiIntegrations ) {
        int fact = sc.param.general.veryVerbose ? 100 : 10;
        int mod = std::max(sc.param.prop.outputCount / fact, 1u);
        if ( callCount % mod == 0 ){
            if ( sc.param.mpi.rank == sc.param.mpi.root ){
                Logger::print("%g fs %u/%u [~%2.1f%%]\n", atomicUnits::to_fs(t), callCount,
                                       sc.param.prop.outputCount, 100 * (double)(callCount) / sc.param.prop.outputCount);
            }
        }
    }
    callCount++;
}

void ObserverBase::reset()
{
    expValues.clear();
    callCount = 0;
}

bool ObserverBase::getStoreRegionDensity(unsigned outputIndex) const
{
    return sc.param.ksr.storeDensityMatrix && (outputIndex % sc.param.ksr.densityOutputStride == 0);
}

ExpectationValuesEntry_t ObserverBase::prepareExpValue(double t, bool storeRegionDensity)
{
    ExpectationValuesEntry_t ev;
    ev.t = t;
    GeomVector3d E = sc.getE(t);
    for(unsigned i=0; i<sc.dim; i++)
        ev.E[i] = E[i];
    for(unsigned i=sc.dim; i<3; i++)
        ev.E[i] = 0;
    ev.region.resize(region.size());
    for(unsigned i=0; i<ev.region.size(); i++){
        ev.region[i].kPointWeightSum = 0;
        for(unsigned dir=0; dir<3; dir++)
            ev.region[i].j[dir] = 0;
        ev.region[i].occupationWan.resize(sc.numWann, 0);
        ev.region[i].occupationHam.resize(sc.numWann, 0);
        if ( storeRegionDensity )
            ev.region[i].wan.resize(sc.numWann * sc.numWann, 0);
    }
    ev.occupationHamMin.resize(sc.numWann, std::numeric_limits<double>::max());
    ev.occupationHamMax.resize(sc.numWann, std::numeric_limits<double>::lowest());
    return ev;
}


ExpectationValues_t ObserverBase::collectExpectationValues() const
{
    ExpectationValues_t rexp;
    rexp.region.resize(regionContext.size());
    for(unsigned i=0; i<regionContext.size(); i++){
        rexp.region[i].name = regionContext[i].name;
        rexp.region[i].dim = regionContext[i].dim;
    }
    for(unsigned ti=0; ti<expValues.size(); ti++){
        const ExpectationValuesEntry_t & ev = expValues[ti];
        rexp.t.push_back(ev.t);
        rexp.E.push_back({ev.E[0], ev.E[1], ev.E[2]});
        rexp.A.push_back({ev.A[0], ev.A[1], ev.A[2]});
        for(unsigned i=0; i<ev.occupationHamMin.size(); i++)
            rexp.occupationHamMin.push_back(ev.occupationHamMin[i]);
        for(unsigned i=0; i<ev.occupationHamMax.size(); i++)
            rexp.occupationHamMax.push_back(ev.occupationHamMax[i]);
        unsigned rid = 0;
        bool storeDensity = false;
        for(unsigned outRegionId=0; outRegionId<regionContext.size(); outRegionId++){
            RegionExpectationValues_t & regionExp = rexp.region[outRegionId];
            const auto &dim = regionContext[outRegionId].dim;
            for(unsigned a=0; a<dim[0]; a++)
                for(unsigned b=0; b<dim[1]; b++)
                    for(unsigned c=0; c<dim[2]; c++){
                        const ExpectationValuesRegion_t &r = ev.region[rid++];
                        regionExp.kPointWeightSum.push_back(r.kPointWeightSum);
                        regionExp.j.push_back({r.j[0], r.j[1], r.j[2]});
                        for(unsigned i=0; i<r.occupationWan.size(); i++)
                            regionExp.occupationWan.push_back(r.occupationWan[i]);
                        for(unsigned i=0; i<r.occupationHam.size(); i++)
                            regionExp.occupationHam.push_back(r.occupationHam[i]);
                        for(unsigned i=0; i<r.wan.size(); i++)
                            regionExp.wan.push_back(r.wan[i]);
                        storeDensity |= r.wan.size() > 0;
                    }
        }
        if ( storeDensity )
            rexp.tDensIndices.push_back(ti);
    }
    return rexp;
}
