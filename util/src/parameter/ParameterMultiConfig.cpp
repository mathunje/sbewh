#include "ParameterMultiConfig.h"

ParameterMultiConfig::ParameterMultiConfig(const std::string name,
                                           const std::string description, std::function<bool(std::string)> callback,
                                           const std::string fixedInputDir) :
    ParameterBase(name, description), callback(callback), fixedInputDir(fixedInputDir)
{
    fname = "";
}

ParameterMultiConfig::~ParameterMultiConfig()
{

}

void ParameterMultiConfig::parseData(const char * str)
{
    fname = std::string(str);
    if ( fixedInputDir.size() > 0 )
        fname = std::filesystem::path(fixedInputDir) / std::filesystem::path(fname).filename();
}

void ParameterMultiConfig::printValue(FILE *f) const
{
    fprintf(f, "%s", fname.c_str());
}

const std::string& ParameterMultiConfig::getFileName() const
{
    return fname;
}
