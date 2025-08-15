#ifndef SBE_WH_OBSERVER_BASE_H
#define SBE_WH_OBSERVER_BASE_H

#include "Parameter.h"
#include "StateContext.h"
#include "ExpectationValues.h"

#include<cmath>
#include<limits>
#include<vector>



struct ExpectationValuesRegion_t {
    double kPointWeightSum;
    double j[3]; // normed by kPointWeightSum
    std::vector<double> occupationHam; // normed by kPointWeightSum
    std::vector<double> occupationWan; // normed by kPointWeightSum
    std::vector< std::complex<double> > wan; // interpolated Wannier density matrix normed by kPointWeightSum
};

struct ExpectationValuesEntry_t {
    double t;
    double E[3];
    double A[3];
    std::vector<ExpectationValuesRegion_t> region;
    std::vector<double> occupationHamMin;
    std::vector<double> occupationHamMax;
};

struct ObserverRegion_t {
    double k[3];
    bool movingFrame;
    double invSigma2;
};

struct RegionContext_t {
    std::string name;
    std::array<unsigned, 3> dim;
};

class ObserverBase {
protected:
    const StateContext &sc;
    unsigned callCount;
    std::vector<ObserverRegion_t> region;
    std::vector<RegionContext_t> regionContext;
    std::vector<ExpectationValuesEntry_t> expValues;

    void basicLog(double t);
    bool getStoreRegionDensity(unsigned outputIndex) const;
    ExpectationValuesEntry_t prepareExpValue(double t, bool storeRegionDensity);

public:
    ObserverBase(const StateContext &sc);
    ExpectationValues_t collectExpectationValues() const;
    void reset();
};

#endif
