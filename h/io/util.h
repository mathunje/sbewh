#ifndef SBE_WH_IO_UTIL_H
#define SBE_WH_IO_UTIL_H


#include<filesystem>

#include "Logger.h"
#include "Parameter.h"


std::string getOutputDirectoryPath(const Parameter_t &param, const std::string fname);
bool createOutputDirectory(Parameter_t &param);
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
