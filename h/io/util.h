#ifndef SBE_WH_IO_UTIL_H
#define SBE_WH_IO_UTIL_H


#include<filesystem>

#include "Logger.h"
#include "Parameter.h"


std::string getOutputDirectoryPath(const Parameter_t &param, const std::string fname);
bool createOutputDirectory(Parameter_t &param);
std::vector<std::string> loadArgsFromInputDirectory(const std::string &dir, bool *ok);

#endif
