#ifndef SBE_WH_IO_UTIL_H
#define SBE_WH_IO_UTIL_H


#include<filesystem>

#include "util/Logger.h"
#include "util/GeomVector.hpp"
#include "OutputParameter.h"


std::string getOutputDirectoryPath(const OutputParameter_t &out, const std::string fname);
bool createOutputDirectory(const OutputParameter_t &out);
std::vector<std::string> loadArgsFromInputDirectory(const std::string &dir, bool *ok);

template<unsigned M, long unsigned N, typename T>
std::array<T, N*M> unravel(const std::array<GeomVector<T, M>, N> &v){
    std::array<T, N*M> res;
    for(unsigned i=0; i<N; i++)
        for(unsigned j=0; j<M; j++)
            res[i*M+j] = v[i].at(j);
    return res;
}

#endif
