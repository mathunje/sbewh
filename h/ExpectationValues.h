#ifndef SBE_WH_EXPECTATION_VALUES_H
#define SBE_WH_EXPECTATION_VALUES_H

#include<array>
#include<vector>
#include<string>


struct RegionExpectationValues_t
{
    std::string name;
    std::array<unsigned, 3> dim;
    std::vector<std::array<double, 3>> j; // dim x t
    std::vector<double> kPointWeightSum; // dim x t
    std::vector<double> occupationHam; // dim x t x occHam
    std::vector<double> occupationWan; // dim x t x occWann
    std::vector< std::complex<double> > wan; // dim x tDens x numWann x numWann
};

struct ExpectationValues_t {
    std::vector<double> t;
    std::vector<unsigned> tDensIndices;
    std::vector<std::array<double, 3>> E;
    std::vector<std::array<double, 3>> A;
    std::vector<double> occupationHamMin; // t x numWann
    std::vector<double> occupationHamMax; // t x numWann
    std::vector<RegionExpectationValues_t> region;
};


#endif
